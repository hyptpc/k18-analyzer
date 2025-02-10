#!/usr/bin/env python3

import argparse
import logging
import logging.config
import numpy as np
import os
import sys
import yaml

import ROOT

import dctdc
import dcdrft
from detector import dc_constants as dcconst
import macrohelper as mh

logger = logging.getLogger(__name__)
ROOT.gStyle.SetHistMinimumZero()

#______________________________________________________________________________
@mh.update_canvas()
def tdctot(c1, name, key):
  hname = f'{name}_{key}_layer0{mh.beamflag_for_param}'
  logger.info(f'hname={hname}')
  h1 = mh.get(hname)
  h1.Draw()
  fig_path = c1.GetTitle().replace('.pdf', f'_{name}_{key}.png')
  ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas()
def tdcfit(c1, name):
  logger.info(f'name={name}')
  hname = f'{name}_CTDC_layer0{mh.beamflag_for_param}'
  h1 = mh.get(hname)
  xmax = h1.GetBinCenter(h1.GetMaximumBin())
  stddev = h1.GetStdDev()
  tdcrange = (xmax - stddev, xmax + stddev)
  params = np.ndarray(4, dtype='float64')
  params[0] = h1.GetMaximum(); params[1] = xmax + 20;
  params[2] = -10; params[3] = 0
  limits = [(h1.GetMaximum()*0.95, h1.GetMaximum()*1.05),
            (xmax, xmax+50),
            (-20, -1),
            (0, h1.GetMaximum()*0.01)]
  h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
  h2 = h1.Clone()
  f1, t0 = dctdc.fit(h1, params=params, limits=limits,
                     drawline=False, drawtext=False)
  m0 = h1.GetMaximum()
  a1 = ROOT.TArrow()
  a1.SetLineColor(ROOT.kRed+1)
  a1.SetFillColor(ROOT.kRed+1)
  a1.DrawArrow(t0, m0*0.11, t0, m0*0.1, 0.04, '|>')
  tex = ROOT.TLatex()
  # tex.SetNDC()
  tex.SetTextAlign(21)
  # tex.SetTextColor(ROOT.kRed+1)
  # tex.SetTextSize(0.07)
  # x = 0.48; y = 0.40
  x = t0; y = m0*0.18
  tex.DrawLatex(x, y, f'p0')
  fig_path = c1.GetTitle().replace('.pdf', f'_{name}_CTDCFit.png')
  ROOT.gPad.Print(fig_path)
  mean = f1.GetParameter(1)
  sigma = abs(f1.GetParameter(2))
  h2.GetXaxis().SetRangeUser(mean-10*sigma, mean+20*sigma)
  h2.Draw()
  fig_path = c1.GetTitle().replace('.pdf', f'_{name}_TDCFitZoom.png')
  ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def drift_time_integral(c1, name):
  logger.info(f'name={name}')
  statx = ROOT.gStyle.GetStatX()
  staty = ROOT.gStyle.GetStatY()
  ROOT.gStyle.SetStatX(0.88)
  ROOT.gStyle.SetStatY(0.6)
  hname = f'{name}_Hit_DriftTime_layer0{mh.beamflag_for_param}'
  h1 = mh.get(hname)
  ret = dcdrft.fit_integral(h1)
  ROOT.gStyle.SetStatX(statx)
  ROOT.gStyle.SetStatY(staty)
  fig_path = c1.GetTitle().replace('.pdf', f'_{name}_DTIntegral.png')
  ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  for n, v in dcconst.items():
    tdctot(n, 'TDC')
    tdctot(n, 'TOT')
    tdctot(n, 'CTDC')
    tdctot(n, 'CTOT')
    tdctot(n, 'HitPat')
    tdctot(n, 'CHitPat')
    tdctot(n, 'Multi')
    tdctot(n, 'CMulti')
    tdctot(n, 'Hit_DriftTime')
    tdctot(n, 'Hit_DriftLength')
    tdctot(n, 'Hit_HitPat')
    tdctot(n, 'Hit_Multi')
    tdcfit(n)
    drift_time_integral(n)
  mh.finalize()

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parser.add_argument('--update', '-u', action='store_true',
                      help='update HodoParam')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(os.path.dirname(__file__), 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  mh.run(parsed.run_list, single_run)

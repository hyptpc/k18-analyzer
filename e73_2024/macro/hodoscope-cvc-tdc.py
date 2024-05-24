#!/usr/bin/env python3

import argparse
import logging
import logging.config
import numpy as np
import os
import sys
import yaml

import ROOT

import hdprm
import macrohelper as mh

logger = logging.getLogger(__name__)
name = 'CVC'
cid = 14
n_seg = 10

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def tdc(ud, tdcrange=(1.21e6, 1.24e6), fit=True):
  logger.info(f'ud={ud}, tdcrange={tdcrange}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(5, 2)
  result_dict = dict()
  for seg in range(n_seg):
    c1.cd(seg+1) #.SetLogy()
    h1 = ROOT.gFile.Get(name + f'_TDC_seg{seg}{ud}{mh.beamflag_for_param}')
    if h1:
      logger.debug(name + f'_TDC_seg{seg}{ud}{mh.beamflag_for_param}')
      h1.RebinX(4)
      if h1.GetEntries() < 1e4:
        h1.RebinX(2)
      # h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
      if fit:
        mean = h1.GetBinCenter(h1.GetMaximumBin())
        sigma = min(4e3, h1.GetStdDev())
        params = np.ndarray(3, dtype='float64')
        params[0] = h1.GetMaximum()
        params[1] = mean
        params[2] = sigma
        limits = [
          (0, h1.GetMaximum()*10),
          (mean - 3*sigma, mean + 3*sigma),
          (1e2, 1e4)
        ]
        result = mh.fit_gaus(h1, params=params, limits=limits)
        key = (cid, 0, seg, 1, 0 if ud == 'U' else 1)
        result_dict[key] = (result.GetParameter(1), -0.000939002)
      else:
        h1.Draw()
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def time(ud='', timerange=(-10, 10), key='Time'):
  logger.info(f'ud={ud}, key={key}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(5, 2)
  pcolor = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kRed+2]
  for seg in range(n_seg):
    c1.cd(seg+1) #.SetLogy()
    for j, b in enumerate(mh.beamflag):
      hname = name + f'_Hit_{key}_seg{seg}{ud}{b}'
      h1 = ROOT.gFile.Get(hname)
      if h1:
        logger.debug(hname)
        h1.SetLineColor(pcolor[j])
        h1.GetXaxis().SetRangeUser(timerange[0], timerange[1])
        h1.Draw('same')
      else:
        logger.warning(f'cannot find {hname}')
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)

#______________________________________________________________________________
def time2d():
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(2, 2)
  keys = ('MeanTime',)
  ranges = ((-10, 10), (-10, 10))
  for i, key in enumerate(keys):
    c1.cd(i+1) #.SetLogy()
    hname = name + f'_Hit_{key}{mh.beamflag_for_param}'
    h1 = ROOT.gFile.Get(hname)
    if h1:
      logger.debug(hname)
      h1.GetXaxis().SetRangeUser(ranges[i][0], ranges[i][1])
      h1.Draw()
    else:
      logger.warning(f'cannot find {hname}')
    c1.cd(i+1+len(keys)).SetLogz()
    hname = name + f'_Hit_{key}_vs_HitPat{mh.beamflag_for_param}'
    h2 = ROOT.gFile.Get(hname)
    if h2:
      logger.debug(hname)
      h2.GetYaxis().SetRangeUser(ranges[i][0], ranges[i][1])
      h2.Draw('colz')
    else:
      logger.warning(f'cannot find {hname}')
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    ret = tdc(ud=ud)
    result_dict.update(ret)
  for key in ['Time',]:
    for ud in ['U', 'D', '']:
      if ud == '':
        key = key.replace('Time', 'MeanTime')
      ret = time(ud=ud, key=key)
  time2d()
  hdprm.output_result(run_info, result_dict, update=parsed.update)
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

#!/usr/bin/env python3

import argparse
import logging
import logging.config
import numpy as np
import os
import sys
import yaml

import ROOT

from detector import hodoscope_constants as hconst
import hdphc
import macrohelper as mh

logger = logging.getLogger(__name__)

#______________________________________________________________________________
@mh.update_canvas()
def t0_adc(c1, hname = f'T0_ADC_seg2U{mh.beamflag_for_param}'):
  logger.info(f'hname={hname}')
  ROOT.gPad.SetLogy()
  h1 = mh.get(hname)
  if h1:
    h1.GetXaxis().SetRangeUser(0, 1000)
    h1.Draw()
    fig_path = c1.GetTitle().replace('.pdf', '_T0_ADC.png')
    ROOT.gPad.Print(fig_path)
    a1 = ROOT.TArrow()
    a1.SetLineColor(ROOT.kRed+1)
    a1.SetFillColor(ROOT.kRed+1)
    p0 = h1.GetBinCenter(h1.GetMaximumBin())
    m0 = h1.GetMaximum()
    h1.GetXaxis().SetRangeUser(p0+40, 1000)
    p1 = h1.GetBinCenter(h1.GetMaximumBin())
    m1 = h1.GetMaximum()
    h1.GetXaxis().SetRangeUser(0, 1000)
    a1.DrawArrow(p0, m0*1.2, p0, m0*1.1, 0.04, '|>')
    a1.DrawArrow(p1, m1*1.2, p1, m1*1.1, 0.04, '|>')
    t1 = ROOT.TLatex()
    t1.SetTextAlign(21)
    t1.DrawLatex(p0, m0*3, 'p0')
    t1.DrawLatex(p1, m1*3, 'p1')
    fig_path = c1.GetTitle().replace('.pdf', '_T0_ADC_HDPRM.png')
    ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas()
def bht_tot(c1):
  ROOT.gPad.SetLogy(0)
  hname = f'BHT_TOT_seg30U{mh.beamflag_for_param}'
  h1 = mh.get(hname)
  h1.Draw()
  fig_path = c1.GetTitle().replace('.pdf', '_BHT_RAW_TOT.png')
  ROOT.gPad.Print(fig_path)
  hname = f'BHT_Hit_TOT_seg30U{mh.beamflag_for_param}'
  h1 = mh.get(hname)
  h1.GetXaxis().SetRangeUser(0, 60)
  h1.Draw()
  fig_path = c1.GetTitle().replace('.pdf', '_BHT_TOT.png')
  ROOT.gPad.Print(fig_path)
  a1 = ROOT.TArrow()
  a1.SetLineColor(ROOT.kRed+1)
  a1.SetFillColor(ROOT.kRed+1)
  p0 = h1.GetBinCenter(h1.GetMaximumBin())
  m0 = h1.GetMaximum()
  a1.DrawArrow(p0, m0*1.05, p0, m0*1.01, 0.04, '|>')
  t1 = ROOT.TLatex()
  t1.SetTextAlign(21)
  t1.DrawLatex(p0, m0*1.08, 'p1')
  fig_path = c1.GetTitle().replace('.pdf', '_BHT_TOT_HDPRM.png')
  ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas()
def hodo_tdc(c1, hname, opt):
  logger.info(f'hname={hname}')
  ROOT.gPad.SetLogy()
  h1 = mh.get(hname)
  if h1:
    h1.Draw()
    fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_TDC.png')
    ROOT.gPad.Print(fig_path)
    a1 = ROOT.TArrow()
    a1.SetLineColor(ROOT.kRed+1)
    a1.SetFillColor(ROOT.kRed+1)
    p0 = h1.GetBinCenter(h1.GetMaximumBin())
    m0 = h1.GetMaximum()
    h1.GetYaxis().SetRangeUser(0.2, m0*8)
    a1.DrawArrow(p0, m0*1.2, p0, m0*1.1, 0.04, '|>')
    t1 = ROOT.TLatex()
    t1.SetTextAlign(21)
    t1.DrawLatex(p0, m0*3, 'p0')
    fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_TDC_HDPRM.png')
    ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas()
def hodo_de(c1, hname, opt):
  logger.info(f'hname={hname}')
  ROOT.gPad.SetLogy(0)
  h1 = mh.get(hname)
  if h1:
    h1.RebinX(4)
    h1.Draw()
    fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_DeltaE.png')
    ROOT.gPad.Print(fig_path)
    a1 = ROOT.TArrow()
    a1.SetLineColor(ROOT.kRed+1)
    a1.SetFillColor(ROOT.kRed+1)
    p0 = h1.GetBinCenter(h1.GetMaximumBin())
    m0 = h1.GetMaximum()
    #h1.GetXaxis().SetRangeUser(0, 1000)
    a1.DrawArrow(p0, m0*1.05, p0, m0*1.01, 0.04, '|>')
    t1 = ROOT.TLatex()
    t1.SetTextAlign(21)
    t1.DrawLatex(p0, m0*1.08, 'MIP')
    fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_DeltaE_HDPRM.png')
    ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas()
def hodo_time(c1, hname, opt, logy=False, label=''):
  logger.info(f'hname={hname}')
  ROOT.gPad.SetLogy(logy)
  h1 = mh.get(hname)
  if h1:
    h1.RebinX(4)
    h1.Draw()
    fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_Time.png')
    ROOT.gPad.Print(fig_path)
    a1 = ROOT.TArrow()
    a1.SetLineColor(ROOT.kRed+1)
    a1.SetFillColor(ROOT.kRed+1)
    p0 = h1.GetBinCenter(h1.GetMaximumBin())
    m0 = h1.GetMaximum()
    #h1.GetXaxis().SetRangeUser(0, 1000)
    h1.GetYaxis().SetRangeUser(0, m0*1.15)
    a1.DrawArrow(p0, m0*1.05, p0, m0*1.01, 0.04, '|>')
    t1 = ROOT.TLatex()
    t1.SetTextAlign(21)
    t1.DrawLatex(p0, m0*1.09, label)
    fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_Time_HDPRM.png')
    ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
@mh.update_canvas()
def hodo_phc(c1, hname, opt):
  logger.info(f'hname={hname}')
  # ROOT.gPad.SetLogz()
  h1 = mh.get(hname)
  h1.Draw('colz')
  fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_Cor.png')
  ROOT.gPad.Print(fig_path)
  params = np.ndarray(3, dtype='float64')
  params[0] = 2
  params[1] = 0
  params[2] = 2
  limits = [
    (0.001, 10),
    (-1, 0.1),
    (-0.5, 10)
  ]
  prof = h1.ProfileX('_pfx', h1.GetYaxis().FindBin(-4),
                     h1.GetYaxis().FindBin(4))
  prof.RebinX(2)
  result = mh.fit_phc(prof, params=params, limits=limits,
                      fitrange=(0.4, 1.8))
  fig_path = c1.GetTitle().replace('.pdf', f'_{opt}_Cor_HDPHC.png')
  ROOT.gPad.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  t0_adc()
  hodo_de(hname=f'T0_Hit_DeltaE_seg2U{mh.beamflag_for_param}', opt='T0')
  hodo_de(hname=f'T0_Hit_DeltaE_seg2{mh.beamflag_for_param}', opt='T0M')
  hodo_tdc(hname=f'T0_TDC_seg2U{mh.beamflag_for_param}', opt='T0')
  hodo_time(hname = f'T0_Hit_Time_seg2U{mh.beamflag_for_param}', opt='T0')
  hodo_time(hname = f'T0_Hit_MeanTime_seg2{mh.beamflag_for_param}', opt='T0M')
  hodo_tdc(hname=f'BHT_TDC_seg30U{mh.beamflag_for_param}', opt='BHT')
  hodo_time(hname=f'BHT_Hit_Time_seg30U{mh.beamflag_for_param}', opt='BHT')
  hodo_time(hname=f'BHT_Hit_MeanTime_seg30{mh.beamflag_for_param}', opt='BHTM')
  bht_tot()
  hodo_de(hname=f'BHT_Hit_DeltaE_seg30U{mh.beamflag_for_param}', opt='BHT')
  hodo_de(hname=f'BHT_Hit_DeltaE_seg30{mh.beamflag_for_param}', opt='BHTM')
  hodo_time(hname=f'T0_seg1_TimeOffset{mh.beamflag_for_param}', opt='T0Ofs', label='p1')
  hodo_phc(hname=f'T0_seg2U_BTOF_vs_DeltaE{mh.beamflag_for_param}', opt='T0')
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

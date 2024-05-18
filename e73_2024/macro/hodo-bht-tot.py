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
import macrohelper

logger = logging.getLogger(__name__)
name = 'BHT'
bht_cid = 1
n_seg = 63
n_seg_one_page = 16
beamflag = '_Pi'

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def tot(start_seg, ud, beamflag='', totrange=(0, 40), fit=True):
  logger.info(f'seg={start_seg}-{start_seg+16}, ud={ud}, beamflag={beamflag}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(4, 4)
  result_dict = dict()
  for i in range(n_seg_one_page):
    c1.cd(i+1) #.SetLogy()
    seg = start_seg + i
    if seg >= n_seg:
      continue
    hname = name + f'_Hit_TOT_seg{seg}{ud}{beamflag}'
    h1 = ROOT.gFile.Get(hname)
    if h1:
      logger.debug(hname)
      h1.GetXaxis().SetRangeUser(totrange[0], totrange[1])
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
          (1, 100)
        ]
        result = macrohelper.fit_gaus(h1, params=params, limits=limits)
        key = (bht_cid, 0, seg, 0, 0 if ud == 'U' else 1)
        result_dict[key] = (0, result.GetParameter(1))
      else:
        h1.Draw()
    else:
      logger.warning(f'cannot find {hname}')
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def de(start_seg, ud='', derange=(0, 3)):
  logger.info(f'seg={start_seg}-{start_seg+16}, ud={ud}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(4, 4)
  pcolor = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kRed+2]
  for i in range(n_seg_one_page):
    c1.cd(i+1) #.SetLogy()
    seg = start_seg + i
    if seg >= n_seg:
      continue
    for j, b in enumerate(['', '_Pi', '_K', '_P']):
      hname = name + f'_Hit_DeltaE_seg{seg}{ud}{b}'
      h1 = ROOT.gFile.Get(hname)
      if h1:
        logger.debug(hname)
        h1.SetLineColor(pcolor[j])
        h1.RebinX(2)
        h1.GetXaxis().SetRangeUser(derange[0], derange[1])
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
  keys = ('MeanTOT', 'DeltaE', )
  ranges = ((0, 30), (0, 3), )
  for i, key in enumerate(keys):
    c1.cd(i+1) #.SetLogy()
    hname = name + f'_Hit_{key}{beamflag}'
    h1 = ROOT.gFile.Get(hname)
    if h1:
      logger.debug(hname)
      h1.GetXaxis().SetRangeUser(ranges[i][0], ranges[i][1])
      h1.Draw()
    else:
      logger.warning(f'cannot find {hname}')
    c1.cd(i+1+len(keys)).SetLogz()
    hname = name + f'_Hit_{key}_vs_HitPat{beamflag}'
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
  macrohelper.initialize(run_info, fig_tail='_bht_tot')
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    for seg in range(4):
      ret = tot(start_seg=seg*16, ud=ud, beamflag=beamflag)
      result_dict.update(ret)
  for ud in ['U', 'D', '']:
    for seg in range(4):
      ret = de(start_seg=seg*16, ud=ud)
  time2d()
  hdprm.output_result(run_info, result_dict, update=parsed.update)
  macrohelper.finalize()

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
  macrohelper.run(parsed.run_list, single_run)

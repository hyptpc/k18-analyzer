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
bht_cid = 1

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def tdc(start_seg, ud, particle='', tdcrange=(1.22e6, 1.26e6), fit=True):
  logger.info(f'seg={start_seg}-{start_seg+16}, ud={ud}, particle={particle}')
  name = 'BHT'
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(4, 4)
  n_seg_one_page = 16
  result_dict = dict()
  for i in range(n_seg_one_page):
    c1.cd(i+1) #.SetLogy()
    seg = start_seg + i
    h1 = ROOT.gFile.Get(name + f'_TDC_seg{seg}{ud}{particle}')
    if h1:
      logger.debug(name + f'_TDC_seg{seg}{ud}{particle}')
      h1.RebinX(4)
      if h1.GetEntries() < 1e4:
        h1.RebinX(2)
      h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
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
        result = macrohelper.fit_gaus(h1, params=params, limits=limits)
        key = (bht_cid, 0, seg, 1, 0 if ud == 'U' else 1)
        result_dict[key] = result
      else:
        h1.Draw()
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  macrohelper.initialize(run_info, fig_tail='_bht')
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    for seg in range(4):
      ret = tdc(start_seg=seg*16, ud=ud, particle='_Pi')
      result_dict.update(ret)
  hdprm.output_result(run_info, result_dict, is_hrtdc=True,
                      update=parsed.update)
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

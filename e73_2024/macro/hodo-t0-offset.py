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
name = 'T0'
cid = 2
n_seg = 5
n_seg_one_page = 16
beamflag = '_Pi'

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def offset(offsetrange=(-2, 2), fit=True):
  logger.info(f'offsetrange={offsetrange}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(3, 2)
  result_dict = dict()
  for seg in range(n_seg):
    c1.cd(seg+1) #.SetLogy()
    hname = name + f'_seg{seg}_TimeOffset{beamflag}'
    h1 = ROOT.gFile.Get(hname)
    if h1:
      logger.debug(hname)
      h1.GetXaxis().SetRangeUser(offsetrange[0], offsetrange[1])
      if fit:
        mean = h1.GetBinCenter(h1.GetMaximumBin())
        sigma = min(0.2, h1.GetStdDev())
        params = np.ndarray(3, dtype='float64')
        params[0] = h1.GetMaximum()
        params[1] = mean
        params[2] = sigma
        limits = [
          (0, h1.GetMaximum()*10),
          (mean - 3*sigma, mean + 3*sigma),
          (0.05, 1.0)
        ]
        result = mh.fit_gaus(h1, params=params, limits=limits)
        key = (cid, 0, seg, 1, 2)
        result_dict[key] = (result.GetParameter(1), f'{1:.6f}')
      else:
        h1.Draw()
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, fig_tail='_t0_offset')
  result_dict = {'generator': os.path.basename(__file__)}
  ret = offset()
  result_dict.update(ret)
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

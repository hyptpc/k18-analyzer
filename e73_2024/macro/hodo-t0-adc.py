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
name = 'T0'
cid = 2
n_seg = 5
n_seg_one_page = 16
beamflag = '_Pi'

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def adc(ud, adcrange=(0, 1000), head='ADC', fit=True):
  logger.info(f'ud={ud}, adcrange={adcrange}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(3, 2)
  result_dict = dict()
  for seg in range(n_seg):
    c1.cd(seg+1) #.SetLogy()
    h1 = ROOT.gFile.Get(name + f'_{head}_seg{seg}{ud}{beamflag}')
    if h1:
      logger.debug(name + f'_{head}_seg{seg}{ud}{beamflag}')
      if h1.GetEntries() < 1e4:
        h1.RebinX(2)
      h1.GetXaxis().SetRangeUser(adcrange[0], adcrange[1])
      if fit:
        mean = h1.GetBinCenter(h1.GetMaximumBin())
        sigma = min(20, h1.GetStdDev())
        params = np.ndarray(3, dtype='float64')
        params[0] = h1.GetMaximum()
        params[1] = mean
        params[2] = sigma
        limits = [
          (0, h1.GetMaximum()*10),
          (mean - 3*sigma, mean + 3*sigma),
          (1, 100)
        ]
        fitrange = (-2, 2) if head == 'AwoT' else (-2, 1)
        result = macrohelper.fit_gaus(h1, params=params, limits=limits,
                                      fitrange=fitrange)
        key = (cid, 0, seg, 0, 0 if ud == 'U' else 1)
        result_dict[key] = result.GetParameter(1)
      else:
        h1.Draw()
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  macrohelper.initialize(run_info, fig_tail='_t0_adc')
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    ret = adc(head='AwoT', ud=ud, adcrange=(0, 250))
    for k, v in ret.items():
      result_dict[k] = [v, 0]
    ret = adc(head='AwT', ud=ud)
    for k, v in ret.items():
      result_dict[k][1] = v
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
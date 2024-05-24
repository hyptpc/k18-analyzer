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
name = 'T1'
cid = 4
n_seg = 1

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def adc(adcrange=(0, 1000), key='ADC', fit=True):
  logger.info(f'adcrange={adcrange}, key={key}, fit={fit}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(2, 2)
  result_dict = dict()
  for iud, ud in enumerate(['U', 'D']):
    for seg in range(n_seg):
      c1.cd(iud*n_seg+seg+1) #.SetLogy()
      flag = '' if key == 'AwoT' else mh.beamflag_for_param
      hname = name + f'_{key}_seg{seg}{ud}{flag}'
      h1 = ROOT.gFile.Get(hname)
      if h1:
        logger.debug(hname)
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
          fitrange = (-2, 2) if key == 'AwoT' else (-2, 1)
          result = mh.fit_gaus(h1, params=params, limits=limits,
                               fitrange=fitrange)
          k = (cid, 0, seg, 0, 0 if ud == 'U' else 1)
          result_dict[k] = result.GetParameter(1)
        else:
          h1.Draw()
      else:
        logger.warning(f'cannot find {hname}')
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def de(derange=(0, 5), key='DeltaE'):
  logger.info(f'derange={derange}, key={key}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(3, 2)
  for iud, ud in enumerate(['U', 'D', '']):
    for seg in range(n_seg):
      c1.cd(iud*n_seg+seg+1) #.SetLogy()
      for j, b in enumerate(mh.beamflag):
        hname = name + f'_Hit_{key}_seg{seg}{ud}{b}'
        h1 = ROOT.gFile.Get(hname)
        if h1:
          logger.debug(hname)
          h1.SetLineColor(mh.beamcolor[j])
          h1.RebinX(5)
          h1.GetXaxis().SetRangeUser(derange[0], derange[1])
          h1.Draw('same')
        else:
          logger.warning(f'cannot find {hname}')
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info=run_info, macro_path=__file__)
  result_dict = {'generator': os.path.basename(__file__)}
  ret = adc(key='AwoT', adcrange=(0, 250))
  for k, v in ret.items():
    result_dict[k] = [v, 0]
  ret = adc(key='AwT')
  for k, v in ret.items():
    result_dict[k][1] = v
  de()
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

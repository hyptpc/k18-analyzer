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
import hdprm
import macrohelper as mh

logger = logging.getLogger(__name__)
name = 'T0'
nseg = hconst[name]['nseg']
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(3, 2))
def adc(c1, ud, adcrange=(0, 1000), key='ADC', fit=True):
  logger.info(f'ud={ud}, adcrange={adcrange}')
  result_dict = dict()
  for seg in range(nseg):
    c1.cd(seg+1) #.SetLogy()
    hname = f'{name}_{key}_seg{seg}{ud}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
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
        k = (hconst[name]['id'], 0, seg, 0, 0 if ud == 'U' else 1)
        result_dict[k] = result.GetParameter(1)
      else:
        h1.Draw()
  return result_dict

#______________________________________________________________________________
@mh.update_canvas(divisions=(3, 2))
def de(c1, ud='', derange=(0, 5), key='DeltaE'):
  logger.info(f'ud={ud}, key={key}')
  pcolor = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kRed+2]
  for seg in range(nseg):
    c1.cd(seg+1) #.SetLogy()
    for j, b in enumerate(mh.beamflag):
      hname = f'{name}_Hit_{key}_seg{seg}{ud}{b}'
      h1 = mh.get(hname)
      if h1:
        h1.SetLineColor(pcolor[j])
        h1.RebinX(5)
        h1.GetXaxis().SetRangeUser(derange[0], derange[1])
        h1.Draw('same')

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    ret = adc(key='AwoT', ud=ud, adcrange=(0, 250))
    for k, v in ret.items():
      result_dict[k] = [v, 0]
    ret = adc(key='AwT', ud=ud)
    for k, v in ret.items():
      result_dict[k][1] = v
  for ud in ['U', 'D', '']:
    de(ud=ud)
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

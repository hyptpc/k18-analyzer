#!/usr/bin/env python3

import argparse
import logging
import logging.config
import numpy as np
import os
import sys
import yaml

import ROOT

from detector import dc_constants as dcconst
import dcdrft
import macrohelper as mh

logger = logging.getLogger(__name__)
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def drift_time(c1, name, key, nlayer=8):
  logger.info(f'name={name}, key={key}, nlayer={nlayer}')
  result_dict = dict()
  for i in range(nlayer):
    c1.cd(i+1) #.SetLogy()
    if key != '':
      ROOT.gPad.SetLogz()
    hname = f'{name}_Hit_DriftTime{key}_layer{i}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
      h1.Draw('colz')
  return result_dict

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def drift_time_integral(c1, name, cid, nlayer=8):
  logger.info(f'name={name}, cid={cid}, nlayer={nlayer}')
  statx = ROOT.gStyle.GetStatX()
  staty = ROOT.gStyle.GetStatY()
  ROOT.gStyle.SetStatX(0.88)
  ROOT.gStyle.SetStatY(0.6)
  result_dict = dict()
  for i in range(nlayer):
    c1.cd(i+1) #.SetLogy()
    hname = f'{name}_Hit_DriftTime_layer{i}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
      ret = dcdrft.fit_integral(h1)
      key = (cid, i, 0)
      result_dict[key] = ret
  ROOT.gStyle.SetStatX(statx)
  ROOT.gStyle.SetStatY(staty)
  return result_dict

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def drift_length(c1, name, key, nlayer=8):
  logger.info(f'name={name}, key={key}, nlayer={nlayer}')
  result_dict = dict()
  for i in range(nlayer):
    c1.cd(i+1) #.SetLogy()
    if len(key) > 0:
      ROOT.gPad.SetLogz()
    hname = f'{name}_Hit_DriftLength{key}_layer{i}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
      h1.Draw('colz')
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for n, v in dcconst.items():
    for key in ['', '_vs_HitPat']:
      drift_time(n, key=key)
    result_dict.update(drift_time_integral(n, v['id']))
    for key in ['', '_vs_HitPat']:
      drift_length(n, key=key)
  dcdrft.output_result(run_info, result_dict, parsed.update)
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

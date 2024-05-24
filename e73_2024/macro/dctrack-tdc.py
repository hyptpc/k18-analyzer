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
import macrohelper as mh

logger = logging.getLogger(__name__)

#______________________________________________________________________________
def tdc(name, nlayer=0, tdcrange=None):
  logger.info(f'name={name}, nlayer={nlayer}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(4, 2)
  for i in range(nlayer):
    c1.cd(i+1) #.SetLogy()
    h1 = ROOT.gFile.Get(f'{name}_CTDC_layer{i}{mh.beamflag_for_param}')
    if h1:
      if tdcrange is None:
        xmax = h1.GetBinCenter(h1.GetMaximumBin())
        stddev = h1.GetStdDev()
        tdcrange = (xmax - stddev, xmax + stddev)
      params = np.ndarray(4, dtype='float64')
      params[0] = -h1.GetMaximum(); params[1] = xmax + 20;
      params[2] = 10; params[3] = 0
      h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
      dctdc.fit(h1, params=params)
      h1.Draw()
  c1.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  if os.path.basename(run_info['bin']) != 'DCTracking':
    logger.error(f'bin must be DCTracking: run_info={run_info}')
    return
  mh.initialize(run_info, fig_tail='_tdc')
  tdc('BLC1a', nlayer=8)
  tdc('BLC1b', nlayer=8)
  tdc('BLC2a', nlayer=8)
  tdc('BLC2b', nlayer=8)
  tdc('BPC1', nlayer=8)
  tdc('BPC2', nlayer=8)
  mh.finalize()

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(os.path.dirname(__file__), 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  mh.run(parsed.run_list, single_run)

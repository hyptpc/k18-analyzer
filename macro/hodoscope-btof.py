#!/usr/bin/env python3

import argparse
import logging
import logging.config
import numpy as np
import os
import sys
import yaml

import ROOT

import hdphc
import macrohelper as mh

logger = logging.getLogger(__name__)
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(3, 2))
def btof_cluster(c1, pair):
  for i, key in enumerate(pair):
    for j, b in enumerate(mh.beamflag):
      c1.cd(i*3+j+1).SetLogz()
      hname = f'{key}{b}'
      h1 = mh.get(hname)
      if h1:
        if 'vs' in key:
          h1.Draw('colz')
        else:
          h1.Draw('colz')
          if j > 0 and h1.GetEntries() > 1e3:
            params = np.ndarray(3, dtype='float64')
            params[0] = h1.GetMaximum()
            params[1] = h1.GetBinCenter(h1.GetMaximumBin())
            params[2] = 0.2
            mh.fit_gaus(h1, params=params, autozoom=False)

#______________________________________________________________________________
@mh.update_canvas(divisions=(3, 2))
def btof_hit(c1, name):
  for i, key in enumerate(['BTOF', 'CBTOF']):
    for j, b in enumerate(mh.beamflag):
      c1.cd(i*3+j+1).SetLogz()
      hname = f'{name}_{key}_vs_DeltaE{b}'
      h1 = mh.get(hname)
      if h1:
        h1.Draw('colz')

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  for name in ['BHT', 'T0']:
    btof_hit(name)
  for pair in [['CTime0', 'CBtof0'],
               ['CBtof0_vs_deT0Seg', 'CBtof0_vs_deBtof0Seg']]:
    btof_cluster(pair)
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

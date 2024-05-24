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
def ftof_cluster():
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  for pair in [['CTime0', 'CFtof0'],
               ['CFtof0_vs_deT0Seg', 'CFtof0_vs_deFtof0Seg']]:
    c1.Clear()
    c1.Divide(3, 2)
    for i, key in enumerate(pair):
      for j, beamflag in enumerate(mh.beamflag):
        c1.cd(i*3+j+1).SetLogz()
        hname = f'{key}{beamflag}'
        h1 = ROOT.gFile.Get(hname)
        if h1:
          logger.debug(hname)
          if 'vs' in key:
            h1.Draw('colz')
          else:
            h1.Draw('colz')
            if j > 0 and h1.GetEntries() > 1e3:
              params = np.ndarray(3, dtype='float64')
              params[0] = h1.GetMaximum()
              params[1] = h1.GetBinCenter(h1.GetMaximumBin())
              params[2] = 0.2
              mh.fit_gaus(h1, params=params)
        else:
          logger.warning(f'cannot find {hname}')
    c1.Modified()
    c1.Update()
    c1.Print(fig_path)

#______________________________________________________________________________
def ftof_hit():
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  for name in ['CVC', 'T0']:
    c1.Clear()
    c1.Divide(3, 2)
    for i, key in enumerate(['FTOF', 'CFTOF']):
      for j, beamflag in enumerate(mh.beamflag):
        c1.cd(i*3+j+1).SetLogz()
        hname = name + f'_{key}_vs_DeltaE{beamflag}'
        h1 = ROOT.gFile.Get(hname)
        if h1:
          logger.debug(hname)
          h1.Draw('colz')
        else:
          logger.warning(f'cannot find {hname}')
    c1.Modified()
    c1.Update()
    c1.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, fig_tail='_ftof')
  ftof_hit()
  ftof_cluster()
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

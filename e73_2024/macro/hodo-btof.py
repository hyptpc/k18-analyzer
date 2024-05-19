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
import macrohelper

logger = logging.getLogger(__name__)

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def btof_cluster():
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  for pair in [['CTime0', 'CBtof0'],
               ['CBtof0_vs_deT0Seg', 'CBtof0_vs_deBtof0Seg']]:
    c1.Clear()
    c1.Divide(4, 2)
    for i, key in enumerate(pair):
      for j, beamflag in enumerate(['', '_Pi', '_K', '_P']):
        c1.cd(i*4+j+1).SetLogz()
        hname = f'{key}{beamflag}'
        h1 = ROOT.gFile.Get(hname)
        if h1:
          logger.debug(hname)
          if 'vs' in key:
            h1.Draw('colz')
          else:
            h1.Draw('colz')
            if h1.GetEntries() > 1e3:
              macrohelper.fit_gaus(h1)
        else:
          logger.warning(f'cannot find {hname}')
    c1.Modified()
    c1.Update()
    c1.Print(fig_path)

#______________________________________________________________________________
def btof_hit():
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  for name in ['BHT', 'T0']:
    c1.Clear()
    c1.Divide(4, 2)
    for i, key in enumerate(['BTOF', 'CBTOF']):
      for j, beamflag in enumerate(['', '_Pi', '_K', '_P']):
        c1.cd(i*4+j+1).SetLogz()
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
  macrohelper.initialize(run_info, fig_tail='_btof')
  btof_hit()
  btof_cluster()
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

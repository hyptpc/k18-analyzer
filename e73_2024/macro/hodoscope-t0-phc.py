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
name = 'T0'
cid = 2
n_seg = 5

ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
def phc(ud, key, fit=True):
  logger.info(f'ud={ud}, key={key}, beamflag={mh.beamflag_for_param}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(3, 2)
  result_dict = dict()
  for seg in range(n_seg):
    c1.cd(seg+1) #.SetLogy()
    hname = name + f'_seg{seg}{ud}_{key}_vs_DeltaE{mh.beamflag_for_param}'
    h1 = ROOT.gFile.Get(hname)
    if h1:
      logger.debug(name + f'_PHC_seg{seg}{ud}{mh.beamflag_for_param}')
      if fit:
        params = np.ndarray(3, dtype='float64')
        params[0] = 2
        params[1] = 0
        params[2] = 2
        limits = [
          (0.001, 10),
          (-1, 0.1),
          (-0.5, 10)
        ]
        h1.Draw('colz')
        prof = h1.ProfileX()
        prof.RebinX(4)
        result = mh.fit_phc(prof, params=params, limits=limits,
                            fitrange=(0.4, 1.6))
        k = (cid, 0, seg, 0 if ud == 'U' else 1)
        result_dict[k] = (1, 3, result.GetParameter(0),
                          result.GetParameter(1), result.GetParameter(2))
        logger.debug(result_dict[k])
      else:
        h1.Draw('colz')
    else:
      logger.warning(f'cannot find {hname}')
  c1.Modified()
  c1.Update()
  c1.Print(fig_path)
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    ret = phc(ud=ud, key='BTOF')
    result_dict.update(ret)
    phc(ud=ud, key='CBTOF', fit=False)
  hdphc.output_result(run_info, result_dict, update=parsed.update)
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

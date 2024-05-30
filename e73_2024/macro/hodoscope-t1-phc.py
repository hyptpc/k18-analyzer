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
import hdphc
import macrohelper as mh

logger = logging.getLogger(__name__)
name = 'T1'
nseg = hconst[name]['nseg']
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(2, 2))
def phc(key, fit=True, pfrange=(-4, 4), fitrange=(0.6, 2.0)):
  logger.info(f'key={key}, beamflag={mh.beamflag_for_param}, '
              +f'pfrange={pfrange}')
  result_dict = dict()
  for iud, ud in enumerate(['U', 'D']):
    for seg in range(nseg):
      c1.cd(iud*nseg+seg+1) #.SetLogy()
      hname = f'{name}_seg{seg}{ud}_{key}_vs_DeltaE{mh.beamflag_for_param}'
      h1 = mh.get(hname)
      if h1:
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
          prof = h1.ProfileX('_pfx', h1.GetYaxis().FindBin(pfrange[0]),
                             h1.GetYaxis().FindBin(pfrange[1]))
          prof.RebinX(4)
          result = mh.fit_phc(prof, params=params, limits=limits,
                              fitrange=fitrange)
          k = (hconst[name]['id'], 0, seg, 0 if ud == 'U' else 1)
          result_dict[k] = (1, 3, result.GetParameter(0),
                            result.GetParameter(1), result.GetParameter(2))
        else:
          h1.Draw('colz')
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  ret = phc(key='BTOF')
  result_dict.update(ret)
  phc(key='CBTOF', fit=False)
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

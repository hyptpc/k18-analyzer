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
import dctdc
import macrohelper as mh

logger = logging.getLogger(__name__)
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def tdc(c1, name, cid, nplane=8, nwire=32, tdcrange=None):
  logger.info(f'name={name}, nplane={nplane}, nwire={nwire}')
  result_dict = dict()
  for i in range(nplane):
    c1.cd(i+1) #.SetLogy()
    hname = f'{name}_CTDC_plane{i}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
      if tdcrange is None:
        xmax = h1.GetBinCenter(h1.GetMaximumBin())
        stddev = h1.GetStdDev()
        tdcrange = (xmax - stddev, xmax + stddev)
      params = np.ndarray(4, dtype='float64')
      params[0] = h1.GetMaximum(); params[1] = xmax + 20;
      params[2] = -10; params[3] = 0
      limits = [(h1.GetMaximum()*0.95, h1.GetMaximum()*1.05),
                (xmax, xmax+50),
                (-20, -1),
                (0, h1.GetMaximum()*0.1)]
      h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
      f1, t0 = dctdc.fit(h1, params=params, limits=limits)
      for wid in range(nwire):
        key = (cid, i, wid)
        result_dict[key] = (t0, -0.8333333333)
  return result_dict

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def drift_time(c1, name, key, nplane=8):
  logger.info(f'name={name}, nplane={nplane}, key={key}')
  result_dict = dict()
  for i in range(nplane):
    c1.cd(i+1) #.SetLogy()
    if len(key) > 0:
      ROOT.gPad.SetLogz()
    hname = f'{name}_Hit_DriftTime{key}_plane{i}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
      h1.Draw('colz')
  return result_dict

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for n, v in dcconst['BcOut'].items():
    result_dict.update(tdc(n, cid=v['id'], nplane=v['nplane'],
                           nwire=v['nwire']))
    for key in ['', '_vs_HitPat']:
      drift_time(n, key=key)
  dctdc.output_result(run_info, result_dict, parsed.update)
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

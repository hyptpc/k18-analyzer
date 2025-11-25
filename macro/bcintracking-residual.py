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
import dcgeo
import macrohelper as mh

logger = logging.getLogger(__name__)
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def residual(c1, name, key, nplane=8, fit=True):
    if not hasattr(residual, 'ref_off'):
        residual.ref_off = None
    logger.info(f'name={name}, key={key}, nplane={nplane}')
    result_dict = dict()

    sigmas = []        
    temps  = []        

    plane_suffix_all = ['U1','UP1','V1','VP1','U2','UP2','V2','VP2']

    for i in range(nplane):
        c1.cd(i+1)
        hname = f'{name}_Track_Residual{key}_plane{i}{mh.beamflag_for_param}'
        h1 = mh.get(hname)

        if not fit:
            h1.Draw('colz')
            continue

        mean  = h1.GetBinCenter(h1.GetMaximumBin())
        sigma = h1.GetRMS()
        params = [h1.GetMaximum(), mean, sigma]
        limits = [(0, h1.GetMaximum()*10),
                  (mean-3*sigma, mean+3*sigma),
                  (0.01, 0.5)]

        result = mh.fit_gaus(h1, params=params, limits=limits, autozoom=False)
        mu     = result.GetParameter(1)
        sig    = result.GetParameter(2)

        sigmas.append(sig)

        plane_suffix = plane_suffix_all[i]
        namekey = f'{name}-{plane_suffix}'

        if namekey == 'BLC1a-U1':
            residual.ref_off = mu
            mu_value = 0.0
        else:
            mu_value = -1.0 * (mu - residual.ref_off)

        temps.append((namekey, mu_value))

    avg_sigma = float(np.mean(sigmas)) if sigmas else 0.2
    for namekey, mu_value in temps:
        result_dict[namekey] = (mu_value, avg_sigma)

    return result_dict

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for n, v in dcconst['BcIn'].items():
    for key in ['', '_vs_DriftLength']:
      if key == '':
        sub_result = residual(n, key=key,fit=True)
        result_dict.update(sub_result)
      else:
        sub_result = residual(n, key=key,fit=False)
        result_dict.update(sub_result)
  dcgeo.output_result(run_info, result_dict, parsed.update)
  mh.finalize()

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parser.add_argument('--update', '-u', action='store_true',
                      help='update DCGeomParam')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(os.path.dirname(__file__), 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  mh.run(parsed.run_list, single_run)

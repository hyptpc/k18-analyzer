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
name = 'BHT'
nseg = hconst[name]['nseg']
nseg_one_page = 16
ROOT.gStyle.SetOptFit(1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 4))
def tdc(c1, start_seg, ud, beamflag='', tdcrange=(1.22e6, 1.26e6), fit=True):
  logger.info(f'seg={start_seg}-{start_seg+16}, ud={ud}, beamflag={beamflag}')
  result_dict = dict()
  for i in range(nseg_one_page):
    c1.cd(i+1) #.SetLogy()
    seg = start_seg + i
    if seg >= nseg:
      continue
    hname = f'{name}_TDC_seg{seg}{ud}{beamflag}'
    h1 = mh.get(hname)
    if h1:
      h1.RebinX(4)
      if h1.GetEntries() < 1e4:
        h1.RebinX(5)
      if h1.GetEntries() < 1e3:
        h1.RebinX(3)
      h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
      if fit:
        mean = h1.GetBinCenter(h1.GetMaximumBin())
        sigma = min(4e3, h1.GetStdDev())
        params = np.ndarray(3, dtype='float64')
        params[0] = h1.GetMaximum()
        params[1] = mean
        params[2] = sigma
        limits = [
          (0, h1.GetMaximum()*10),
          (mean - 3*sigma, mean + 3*sigma),
          (1e2, 1e4)
        ]
        result = mh.fit_gaus(h1, params=params, limits=limits)
        key = (hconst[name]['id'], 0, seg, 1, 0 if ud == 'U' else 1)
        result_dict[key] = (result.GetParameter(1), -0.0009390020)
      else:
        h1.Draw()
  return result_dict

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 4))
def time(c1, start_seg, ud='', timerange=(-20, 20), key='Time'):
  logger.info(f'seg={start_seg}-{start_seg+16}, ud={ud}, key={key}')
  for i in range(nseg_one_page):
    c1.cd(i+1) #.SetLogy()
    seg = start_seg + i
    if seg >= nseg:
      continue
    for j, b in enumerate(mh.beamflag):
      hname = name + f'_Hit_{key}_seg{seg}{ud}{b}'
      h1 = mh.get(hname)
      if h1:
        h1.SetLineColor(mh.beamcolor[j])
        h1.GetXaxis().SetRangeUser(timerange[0], timerange[1])
        h1.Draw('same')

#______________________________________________________________________________
@mh.update_canvas(divisions=(2, 2))
def time2d(c1):
  keys = ('MeanTime', #'CMeanTime', 'MeanTOT'
          )
  ranges = ((-10, 10), (-10, 10), (0, 30))
  for i, key in enumerate(keys):
    c1.cd(i+1) #.SetLogy()
    hname = name + f'_Hit_{key}{mh.beamflag_for_param}'
    h1 = mh.get(hname)
    if h1:
      h1.GetXaxis().SetRangeUser(ranges[i][0], ranges[i][1])
      h1.Draw()
    c1.cd(i+1+len(keys)).SetLogz()
    hname = name + f'_Hit_{key}_vs_HitPat{mh.beamflag_for_param}'
    h2 = mh.get(hname)
    if h2:
      h2.GetYaxis().SetRangeUser(ranges[i][0], ranges[i][1])
      h2.Draw('colz')

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  result_dict = {'generator': os.path.basename(__file__)}
  for ud in ['U', 'D']:
    for seg in range(4):
      ret = tdc(start_seg=seg*16, ud=ud, beamflag=mh.beamflag_for_param)
      result_dict.update(ret)
  for key in ['Time', # 'CTime'
              ]:
    for ud in ['U', 'D']:
      for seg in range(4):
        ret = time(start_seg=seg*16, ud=ud, key=key)
  for key in ['MeanTime', #'CMeanTime', 'MeanTOT'
              ]:
    for seg in range(4):
      if key == 'MeanTOF':
        timerange=(0, 20)
      else:
        timerange=(-20, 20)
      ret = time(start_seg=seg*16, timerange=timerange, key=key)
  time2d()
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

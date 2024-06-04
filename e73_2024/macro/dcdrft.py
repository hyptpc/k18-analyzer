#!/usr/bin/env python3

import array
import datetime
import logging
import os
import shutil

import ROOT

import conf
import macrohelper as mh

logger = logging.getLogger('__main__').getChild(__name__)
tmp_dir = os.path.join(os.path.dirname(__file__), 'tmp')
os.makedirs(tmp_dir, exist_ok=True)

MaxDriftTime = {
  'BLC1a': 200,
  'BLC1b': 200,
  'BLC2a': 200,
  'BLC2b': 200,
  'BPC2': 200,
  'BPC1': 200,
}

MaxDriftLength = {
  'BLC1a': 4.0,
  'BLC1b': 4.0,
  'BLC2a': 2.5,
  'BLC2b': 2.5,
  'BPC2': 3.6,
  'BPC1': 3.0,
}

#______________________________________________________________________________
def fit_integral(h1, params=None, limits=None):
  name = h1.GetName().split('_')[0]
  layer = h1.GetName().split('layer')[1].split('_')[0]
  max_dt = MaxDriftTime[name]
  max_dl = MaxDriftLength[name]
  bin1 = h1.FindBin(-1)
  bin2 = h1.FindBin(max_dt)
  entry = h1.Integral(bin1, bin2)
  dt = array.array('d')
  dl = array.array('d')
  for i in range(bin1, bin2+1):
    integral = h1.Integral(bin1, i)
    dt.append(h1.GetBinCenter(i))
    dl.append(max_dl*integral/entry)
  g1 = ROOT.TGraph(bin2-bin1+1, dt, dl)
  g1.GetXaxis().SetTitleOffset(1.1)
  g1.GetXaxis().SetLimits(-max_dt*0.1, max_dt*1.2)
  g1.GetYaxis().SetRangeUser(-max_dl*0.1, max_dl*1.2)
  g1.SetNameTitle(h1.GetName().replace('DriftTime', 'DriftFunction')
                  .replace(mh.beamflag_for_param, ''),
                  f'DriftFunction {name} L{layer}; '
                  +f'Drift Time [ns]; Drift Length [mm]')
  g1.SetMarkerStyle(8)
  g1.SetMarkerSize(0.4)
  g1.Draw('APC')

  # f1 = ROOT.TF1('f1', 'pol5')
  # f1.FixParameter(0, 0)
  # f1.SetParameter(1, 0.1)
  # f1.SetParameter(2, -0.02)
  # f1.SetParameter(3, 0.0002)
  # f1.SetParameter(4, -0.000007)
  # f1.SetParameter(5, 0.0000000004)
  # g1.Fit('f1', 'Q', '', -1, max_dt)
  return g1 #, f1

#______________________________________________________________________________
def output_result(run_info, result_dict, update=False):
  logger.debug(run_info)
  dcdrft_path = conf.get(run_info, 'DCDRFT')
  output_path = os.path.join(
    tmp_dir, f'DCDriftParam_{run_info["key"]:05d}.root')
  ref_path = ROOT.gFile.GetName()
  f = ROOT.gFile.Open(output_path, 'recreate')
  ROOT.TNamed('datetime', str(datetime.datetime.now())).Write()
  ROOT.TNamed('reference', ref_path).Write()
  for k, g in result_dict.items():
    if k != 'generator':
      g.Write()
  f.Close()
  if update:
    dcdrft_dir = os.path.dirname(dcdrft_path)
    logger.info(
      f'update {os.path.join(dcdrft_dir, os.path.basename(output_path))}')
    shutil.copy2(output_path, dcdrft_dir)
    conf.replace(run_info, 'DCDRFT', os.path.basename(output_path))

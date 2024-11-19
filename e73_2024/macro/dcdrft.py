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

# DCGeomLayerId = {
#   'BLC1a': [1, 2, 3, 4, 5, 6, 7, 8],
#   'BLC1b': [9, 10, 11, 12, 13, 14, 15, 16],
#   'BLC2a': [1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008],
#   'BLC2b': [1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016],
#   'BPC1': [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008],
#   'BPC2': [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016],
# }

#______________________________________________________________________________
def fit_integral(h1, params=None, limits=None):
  name = h1.GetName().split('_')[0]
  plane = h1.GetName().split('plane')[1].split('_')[0]
  # name.replace('layer' + str(layer), '')
  # layer = DCGeomLayerId[name][int(plane)]
  # name += 'layer' + str(layer)
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
                  # .replace('layer' + str(plane), 'layer' + str(layer))
                  .replace(mh.beamflag_for_param, ''),
                  f'DriftFunction {name} plane{plane}; '
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
    base_path = os.path.join(dcdrft_dir, os.path.basename(output_path))
    logger.info(f'update {base_path}')
    f = ROOT.gFile.Open(base_path, 'update')
    if not f:
      f = ROOT.gFile.Open(base_path, 'create')
    ROOT.TNamed('datetime', str(datetime.datetime.now())).Write(
      '', ROOT.TObject.kOverwrite)
    ROOT.TNamed('reference', ref_path).Write(
      '', ROOT.TObject.kOverwrite)
    for k, g in result_dict.items():
      if k != 'generator':
        g.Write('', ROOT.TObject.kOverwrite)
    f.Close()

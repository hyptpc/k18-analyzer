#!/usr/bin/env python3

import argparse
import logging
import logging.config
import os
import sys
import yaml

import ROOT

import macrohelper as mh

logger = logging.getLogger(__name__)


@mh.update_canvas()
def draw_one_data(c1, name, key):
  hist1d_base = ["Mean", "Max", "RMS", "LocMax", "Min"]
  extra_baseline = ["p0", "p1", "p2"]

  if key == "_Baseline":
    hist_list = hist1d_base + extra_baseline
    nx, ny = 4, 2
  else:
    hist_list = hist1d_base
    nx, ny = 3, 2

  c1.Divide(nx, ny)

  for ipad, hname in enumerate(hist_list, start=1):
    c1.cd(ipad).SetLogy()
    h1 = mh.get(f"{name}_FADC{key}_{hname}")
    if h1:
      if hname == "Max":
        h1.GetXaxis().SetRangeUser(0,2000)  
      h1.Draw()


#______________________________________________________________________________
@mh.update_canvas(divisions=(2, 2))
def draw_one_data_noise(c1, name):
  hist_noise_list = ["Max", "RMSfront", "RMSmiddle", "Adcdiff"]
  for ipad, hname in enumerate(hist_noise_list, start=1):
    c1.cd(ipad).SetLogy()
    h1 = mh.get(f"{name}_FADC_Noise_{hname}")
    if h1:
      if "RMS" in hname:
        h1.GetXaxis().SetRangeUser(0, 100)
      h1.Draw()
    
#______________________________________________________________________________
@mh.update_canvas(divisions=(3, 2))
def draw_one_data2d(c1, name):
  hist2d_list = ["Baseline", "Before", "After", "Good", "Noise"]
  for ipad, hname in enumerate(hist2d_list, start=1):
    c1.cd(ipad)
    h1 = mh.get(f'{name}_FADC_{hname}')
    if h1:
      h1.Draw()

#______________________________________________________________________________
@mh.update_canvas(divisions=(2, 1))
def draw_one_data2dpoly(c1, name):
  hist2dpoly_list = ["Baseline", "Noise"]
  for ipad, hname in enumerate(hist2dpoly_list, start=1):
    c1.cd(ipad).SetLogz()
    h1 = mh.get(f'{name}_HitPat_{hname}')
    if h1:
      if hname == "Baseline":
        h1.SetMinimum(0)
        h1.SetMaximum(10)
      h1.Draw("colz")
      ROOT.gPad.Update()

      palette = h1.GetListOfFunctions().FindObject("palette")
      if palette:
        palette.Draw()

#______________________________________________________________________________
def draw(name):
  logger.info(f'name={name}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  for key in ["","_Cor","_Baseline"]:
    draw_one_data(name, key)
  draw_one_data_noise(name)
  draw_one_data2d(name)
  draw_one_data2dpoly(name)

#______________________________________________________________________________
def single_run(run_info):
  mh.initialize(run_info, __file__)
  draw('TPC')
  mh.finalize()
    
#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(os.path.dirname(__file__), 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  mh.run(parsed.run_list, single_run)

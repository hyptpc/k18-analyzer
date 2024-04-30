#!/usr/bin/env python3

import argparse
import logging
import logging.config
import os
import sys
import ROOT
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

macro_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(macro_dir)
sys.path.append(os.path.join(
  os.path.dirname(macro_dir), 'runmanager'))
sys.path.append(os.path.join(
  os.path.dirname(macro_dir), 'runmanager', 'module'))

import runlist
# import histtool as h

#______________________________________________________________________________
def draw(name, nseg=0, adcdiv=None, tdcdiv=None, trailingdiv=None,
         totdiv=None, ud=False):
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  hitmulti = ['_OR', '_AND'] if ud else ['']
  ud = ['U', 'D'] if ud else ['']

  if adcdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(adcdiv[0], adcdiv[1])
      for i in range(nseg):
        c1.cd(i+1).SetLogy()
        h1 = ROOT.gFile.Get(name + f'_ADC_seg{i}{s}')
        if h1: h1.Draw()
        h2 = ROOT.gFile.Get(name + f'_AWT_seg{i}{s}')
        if h2: h2.SetLineColor(ROOT.kRed+1)
        if h2: h2.Draw('same')
      c1.Print(fig_path)

  if tdcdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(tdcdiv[0], tdcdiv[1])
      for i in range(nseg):
        c1.cd(i+1)
        h1 = ROOT.gFile.Get(name + f'_TDC_seg{i}{s}')
        if h1: h1.Draw()
      c1.Print(fig_path)

  if trailingdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(trailingdiv[0], trailingdiv[1])
      for i in range(nseg):
        c1.cd(i+1)
        h1 = ROOT.gFile.Get(name + f'_Trailing_seg{i}{s}')
        if h1: h1.Draw()
      c1.Print(fig_path)

  if totdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(totdiv[0], totdiv[1])
      for i in range(nseg):
        c1.cd(i+1)
        h1 = ROOT.gFile.Get(name + f'_TOT_seg{i}{s}')
        if h1: h1.Draw()
      c1.Print(fig_path)

  c1.Clear()
  c1.Divide(2, 2)
  i = 1
  for s in hitmulti:
    c1.cd(i)
    h1 = ROOT.gFile.Get(f'{name}_HitPat{s}')
    if h1: h1.Draw()
    i = i + 1
    c1.cd(i)
    h1 = ROOT.gFile.Get(f'{name}_Multi{s}')
    if h1: h1.Draw()
    i = i + 1
  c1.Print(fig_path)

#______________________________________________________________________________
def run(run_info):
  logger.debug(run_info)
  try:
    root_path = run_info['root']
    fig_path = os.path.basename(root_path).replace('.root', '.pdf')
    fig_path = os.path.join(run_info['fig'], fig_path)
  except KeyError as e:
    logger.error(f'KeyError: {e} not found in {run_info}')
    return
  ROOT.gROOT.SetBatch()
  c1 = ROOT.TCanvas('c1', fig_path, 1200, 800)
  f1 = ROOT.TFile(root_path)
  if not f1.IsOpen():
    logger.error('root file is not open.')
    return
  logger.info(f'root_path = {root_path}')
  logger.info(f'fig_path = {fig_path}')
  c1.Print(fig_path+'[')
  draw('TriggerFlag', nseg=32, tdcdiv=(8, 4))
  draw('BHT', nseg=63, tdcdiv=(8, 8), totdiv=(8, 8), ud=True)
  draw('AC', nseg=1, adcdiv=(2, 2), tdcdiv=(2, 2))
  c1.Print(fig_path+']')

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(macro_dir, 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  runlist_manager = runlist.RunlistManager()
  runlist_manager.set_run_list(parsed.run_list)
  for run_info in runlist_manager.get_run_list():
    run(run_info)

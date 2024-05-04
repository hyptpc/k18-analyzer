#!/usr/bin/env python3

import argparse
import logging
import logging.config
import os
import sys
import ROOT
import pandas as pd
import yaml

# myname = os.path.splitext(os.path.basename(__file__))[0]
logger = logging.getLogger(__name__) # .getChild(myname)

macro_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(macro_dir)
sys.path.append(os.path.join(
  os.path.dirname(macro_dir), 'runmanager'))
sys.path.append(os.path.join(
  os.path.dirname(macro_dir), 'runmanager', 'module'))

import runlist

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(1110)

#______________________________________________________________________________
def daq():
  node_list = ['EB', 'FE_VME', 'FE_HUL', 'FE_VEASIROC']
  logger.info(f'node_list={node_list}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  c1.Clear()
  c1.Divide(2, 2)
  for i, head in enumerate(node_list):
    c1.cd(i+1)
    h1 = ROOT.gFile.Get(f'{head}_DataSize')
    if h1:
      h1.Draw('colz')
  c1.Print(fig_path)

#______________________________________________________________________________
def hodo(name, nseg=0, adcdiv=None, tdcdiv=None, trailingdiv=None,
         totdiv=None, ud=True):
  logger.info(f'name={name}, nseg={nseg}, adcdiv={adcdiv}, '
              +'tdcdiv={tdcdiv}, totdiv={totdiv}, ud={ud}')
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
        if h1:
          h1.Draw()
        h2 = ROOT.gFile.Get(name + f'_AWT_seg{i}{s}')
        if h2:
          h2.SetLineColor(ROOT.kRed+1)
          h2.Draw('same')
      c1.Print(fig_path)

  if tdcdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(tdcdiv[0], tdcdiv[1])
      for i in range(nseg):
        c1.cd(i+1).SetLogy()
        h1 = ROOT.gFile.Get(name + f'_TDC_seg{i}{s}')
        if h1:
          h1.RebinX(4)
          h1.Draw()
      c1.Print(fig_path)

  # if trailingdiv is not None:
  #   for s in ud:
  #     c1.Clear()
  #     c1.Divide(trailingdiv[0], trailingdiv[1])
  #     for i in range(nseg):
  #       c1.cd(i+1)
  #       h1 = ROOT.gFile.Get(name + f'_Trailing_seg{i}{s}')
  #       if h1: h1.Draw()
  #     c1.Print(fig_path)

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
    if h1:
      h1.Draw()
      efficiency(h1)
    i = i + 1
  c1.Print(fig_path)

#______________________________________________________________________________
def dc(name, nlayer=0, tdcdiv=None):
  logger.info(f'name={name}, nlayer={nlayer}, tdcdiv={tdcdiv}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  if tdcdiv is not None:
    for htype in ['TDC', 'Trailing', 'TOT', 'HitPat', 'Multi']:
      c1.Clear()
      c1.Divide(tdcdiv[0], tdcdiv[1])
      for i in range(nlayer):
        c1.cd(i+1) #.SetLogy()
        h1 = ROOT.gFile.Get(f'{name}_{htype}_layer{i}')
        if h1:
          h1.Draw()
          if htype == 'Multi':
            efficiency(h1)
        h2 = ROOT.gFile.Get(f'{name}_C{htype}_layer{i}')
        if h2:
          h2.SetLineColor(ROOT.kRed+1)
          h1.GetYaxis().SetRangeUser(
            0, max(h1.GetMaximum(), h2.GetMaximum())*1.05)
          h2.Draw('same')
          if htype == 'Multi':
            efficiency(h2)
        # if htype == 'TDC' or htype == 'Trailing' or htype == 'TOT':
        #   h2 = ROOT.gFile.Get(f'{name}_{htype}1st_layer{i}')
        #   if h2:
        #     h2.SetLineColor(ROOT.kRed+1)
        #     h2.Draw('same')
      c1.Print(fig_path)

#______________________________________________________________________________
def efficiency(h1):
  nof0 = h1.GetBinContent(1)
  nall = h1.GetEntries()
  eff = 1 - nof0/nall
  tex = ROOT.TLatex()
  tex.SetNDC()
  tex.SetTextAlign(11)
  tex.SetTextColor(h1.GetLineColor())
  tex.SetTextSize(0.08)
  x = 0.45
  y = 0.70 if h1.GetLineColor() == 602 else 0.62
  tex.DrawLatex(x, y, f'eff. {eff:.3f}')
  logger.debug(f'{h1.GetTitle()} eff.={eff:.3f}')

#______________________________________________________________________________
def run(run_list):
  logger.debug(f'set {run_list}')
  if not os.path.isfile(run_list):
    logger.error(f'No such file: {run_list}')
    return
  # logger.info(f'run={yaml.safe_load(open(run_list, "r"))["RUN"].keys()}')
  runlist_manager = runlist.RunlistManager()
  runlist_manager.set_run_list(run_list)
  run_list = runlist_manager.get_run_list()
  logger.debug(f'run_list={run_list}')
  for run_info in runlist_manager.get_run_list():
    single_run(run_info)
  logger.info('done')

#______________________________________________________________________________
def single_run(run_info):
  logger.debug(run_info)
  try:
    root_path = run_info['root']
    fig_path = os.path.basename(root_path).replace('.root', '.pdf')
    fig_path = os.path.join(run_info['fig'], fig_path)
  except KeyError as e:
    logger.error(f'KeyError: {e} not found in {run_info}')
    return
  c1 = ROOT.TCanvas('c1', fig_path, 1200, 800)
  f1 = ROOT.TFile(root_path)
  if not f1.IsOpen():
    logger.error('root file is not open.')
    return
  logger.info(f'open {root_path}')
  c1.Print(fig_path+'[')
  hodo('TriggerFlag', nseg=32, tdcdiv=(8, 4), ud=False)
  hodo('BHT', nseg=63, tdcdiv=(8, 8), totdiv=(8, 8))
  hodo('AC', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2), ud=False)
  hodo('T1', nseg=1, adcdiv=(2, 2), tdcdiv=(2, 2))
  hodo('T0', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2))
  hodo('T0new', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2))
  hodo('DEF', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2))
  hodo('Veto', nseg=4, adcdiv=(2, 2), tdcdiv=(2, 2))
  hodo('BTC', nseg=4, adcdiv=(2, 2), tdcdiv=(2, 2))
  # hodo('CDH', nseg=36, adcdiv=(6, 6), tdcdiv=(6, 6))
  # hodo('PbG', nseg=40, adcdiv=(7, 6), tdcdiv=(7, 6))
  # hodo('PbF2', nseg=40, adcdiv=(7, 6), tdcdiv=(7, 6))
  hodo('CVC', nseg=9, adcdiv=(3, 3), tdcdiv=(3, 3))
  hodo('NC', nseg=6, adcdiv=(3, 2), tdcdiv=(3, 2))
  dc('BLC1a', nlayer=8, tdcdiv=(4, 2))
  dc('BLC1b', nlayer=8, tdcdiv=(4, 2))
  dc('BLC2a', nlayer=8, tdcdiv=(4, 2))
  dc('BLC2b', nlayer=8, tdcdiv=(4, 2))
  dc('BPC1', nlayer=8, tdcdiv=(4, 2))
  dc('BPC2', nlayer=8, tdcdiv=(4, 2))
  dc('VFT', nlayer=14, tdcdiv=(5, 3))
  daq()
  c1.Print(fig_path+']')
  logger.info(f'save {fig_path}')

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(macro_dir, 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  run(parsed.run_list)

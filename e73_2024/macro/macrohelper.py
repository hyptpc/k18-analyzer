#!/usr/bin/env python3

import logging
import multiprocessing as mp
import os
import sys

import ROOT

logger = logging.getLogger('__main__').getChild(__name__)

macro_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(macro_dir)
sys.path.append(os.path.join(
  os.path.dirname(macro_dir), 'runmanager'))
sys.path.append(os.path.join(
  os.path.dirname(macro_dir), 'runmanager', 'module'))

import runlist

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(1110)

c1 = None

#______________________________________________________________________________
def daq():
  node_list = ['EB', 'FE_VME', 'FE_HUL', 'FE_VEASIROC']
  logger.info(f'node_list={node_list}')
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
def efficiency(h1):
  nof0 = h1.GetBinContent(1)
  nall = h1.GetEntries()
  if nall == 0:
    logger.warning(f'{h1.GetTitle()} nentries=0')
    return
  eff = 1 - nof0/nall
  tex = ROOT.TLatex()
  tex.SetNDC()
  tex.SetTextAlign(11)
  tex.SetTextColor(h1.GetLineColor())
  tex.SetTextSize(0.08)
  x = 0.45
  if h1.GetLineColor() == ROOT.kBlack:
    y = 0.78
  elif h1.GetLineColor() == ROOT.kBlue+2:
    y = 0.70
  elif h1.GetLineColor() == ROOT.kGreen+2:
    y = 0.62
  elif h1.GetLineColor() == ROOT.kRed+2:
    y = 0.54
  elif h1.GetLineColor() == ROOT.kRed+1:
    y = 0.62
  else:
    y = 0.70
  tex.DrawLatex(x, y, f'eff. {eff:.3f}')
  logger.debug(f'{h1.GetTitle()} eff.={eff:.3f}')

#______________________________________________________________________________
def finalize():
  fig_path = c1.GetTitle()
  c1.Print(fig_path+']')
  logger.info(f'save {fig_path}')

#______________________________________________________________________________
def fit_gaus(h1, params, limits=None):
  f1 = ROOT.TF1('f1', 'gaus')
  f1.SetParameters(params)
  if limits is not None:
    for i, l in enumerate(limits):
      f1.SetParLimits(i, l[0], l[1])
  for j in range(3):
    mean = f1.GetParameter(1)
    sigma = f1.GetParameter(2)
    h1.Fit('f1', 'Q', '', mean - 2*sigma, mean + 2*sigma)
  return f1

#______________________________________________________________________________
def initialize(run_info, fig_tail=''):
  logger.debug(run_info)
  try:
    root_path = run_info['root']
    fig_path = os.path.basename(root_path).replace('.root', fig_tail+'.pdf')
    fig_path = os.path.join(run_info['fig'], fig_path)
  except KeyError as e:
    logger.error(f'KeyError: {e} not found in {run_info}')
    return
  global c1
  c1 = ROOT.TCanvas('c1', fig_path, 1200, 800)
  ROOT.gFile = ROOT.TFile.Open(root_path)
  if ROOT.gFile == None or not ROOT.gFile.IsOpen():
    logger.error('root file is not open.')
    return
  logger.info(f'open {root_path}')
  c1.Print(fig_path+'[')
  title()
  status()

#______________________________________________________________________________
def run(run_list, target):
  logger.debug(f'set {run_list}')
  if not os.path.isfile(run_list):
    logger.error(f'No such file: {run_list}')
    return
  runlist_manager = runlist.RunlistManager()
  runlist_manager.set_run_list(run_list)
  run_list = runlist_manager.get_run_list()
  logger.debug(f'run_list={run_list}')
  proc_list = list()
  for run_info in runlist_manager.get_run_list():
    proc = mp.Process(target=target, args=(run_info,))
    proc.start()
    proc_list.append(proc)
  try:
    for proc in proc_list:
      proc.join()
    logger.info('done')
  except KeyboardInterrupt:
    print()
    for proc in proc_list:
      proc.terminate()
    logger.info('terminated')

#______________________________________________________________________________
def status():
  fig_path = c1.GetTitle()
  h1 = ROOT.gFile.Get('Status')
  if h1:
    prev_opt = ROOT.gStyle.GetOptStat()
    ROOT.gStyle.SetOptStat(10)
    entry = h1.GetBinContent(1)
    passed = h1.GetBinContent(h1.GetNbinsX())
    eff = passed/entry if entry > 0 else ROOT.TMath.QuietNaN()
    logger.info(f'entry={entry:.0f}, passed={passed:.0f} '
                + f'({eff:.2f})')
    c1.Clear()
    h1.Draw()
    c1.Print(fig_path)
    ROOT.gStyle.SetOptStat(prev_opt)

#______________________________________________________________________________
def title():
  fig_path = c1.GetTitle()
  name = ROOT.gFile.GetName()
  now = ROOT.TTimeStamp()
  now.Add(-ROOT.TTimeStamp.GetZoneOffset())
  tex = ROOT.TLatex()
  tex.SetNDC()
  tex.SetTextAlign(22)
  c1.Clear()
  tex.SetTextSize(0.04)
  tex.DrawLatex(0.5, 0.64, str(sys.argv))
  tex.DrawLatex(0.5, 0.58, os.path.dirname(name))
  tex.SetTextSize(0.08)
  tex.DrawLatex(0.5, 0.50, os.path.basename(name))
  tex.DrawLatex(0.5, 0.42, 'Creation time: '+now.AsString('s'))
  c1.Print(fig_path)
  logger.info(sys.argv)
  logger.info(f'time={now.AsString("s")}')
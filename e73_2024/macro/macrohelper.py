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
beamflag = ['', '_Pi', '_K']
beamcolor = [ROOT.kBlack, ROOT.kRed+2, ROOT.kBlue+2]
beamflag_for_param = beamflag[1]

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
  if h1.GetLineColor() == beamcolor[0]:
    y = 0.78
  elif h1.GetLineColor() == beamcolor[1]:
    y = 0.70
  elif h1.GetLineColor() == beamcolor[2]:
    y = 0.62
  elif h1.GetLineColor() == ROOT.kRed+1:
    y = 0.54
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
def fit_gaus(h1, params=None, limits=None, fitrange=(-2, 2), autozoom=True):
  f1 = ROOT.TF1('f1', 'gaus')
  f1.SetLineWidth(1)
  if params is not None:
    f1.SetParameters(params)
  if limits is not None:
    for i, l in enumerate(limits):
      f1.SetParLimits(i, l[0], l[1])
  for j in range(3):
    mean = f1.GetParameter(1)
    sigma = f1.GetParameter(2)
    h1.Fit('f1', 'Q', '', mean + fitrange[0]*sigma, mean + fitrange[1]*sigma)
  if autozoom:
    h1.GetXaxis().SetRangeUser(mean-10*sigma, mean+20*sigma)
  return f1

#______________________________________________________________________________
def fit_phc(h1, params, limits=None, fitrange=(0.5, 1.5)):
  f1 = ROOT.TF1('f1', '-[0]/sqrt(x-[1])+[2]', -1, 10);
  f1.SetLineWidth(1)
  f1.SetParameters(params)
  if limits is not None:
    for i, l in enumerate(limits):
      f1.SetParLimits(i, l[0], l[1])
  h1.Fit('f1', 'Q', 'same', fitrange[0], fitrange[1])
  return f1

#______________________________________________________________________________
def initialize(run_info, macro_path, comment=''):
  logger.debug(run_info)
  try:
    root_path = run_info['root']
    target_bin = os.path.basename(run_info['bin'])
    macro_path = os.path.basename(macro_path).lower()
    if target_bin.lower() not in macro_path:
      logger.error(f'bin must be {target_bin}: run_info={run_info}, '
                   +f'macro_path={macro_path}')
      sys.exit(1)
    fig_tail = (macro_path.replace(target_bin.lower(), '')
                .replace('.py', '').replace('-', '_'))
    fig_path = os.path.basename(root_path).replace('.root', fig_tail+'.pdf')
    fig_path = os.path.join(run_info['fig'], fig_path)
  except KeyError as e:
    logger.error(f'KeyError: {e} not found in {run_info}')
    sys.exit(1)
  global c1
  c1 = ROOT.TCanvas('c1', fig_path, 1200, 800)
  ROOT.gFile = ROOT.TFile.Open(root_path)
  if not ROOT.gFile or not ROOT.gFile.IsOpen():
    logger.error('root file is not open.')
    sys.exit(1)
  logger.info(f'open {root_path}')
  c1.Print(fig_path+'[')
  title(comment=comment)
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
    proc = mp.Process(target=target, args=(run_info,),
                      name=f'[{run_info["key"]:05d}]')
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
    ROOT.gStyle.SetOptStat(0)
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
def title(comment):
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
  tex.SetTextSize(0.04)
  tex.DrawLatex(0.5, 0.34, comment)
  c1.Print(fig_path)
  logger.info(sys.argv)
  logger.info(f'time={now.AsString("s")}')

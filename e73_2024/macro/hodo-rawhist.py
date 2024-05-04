#!/usr/bin/env python3

import argparse
import logging
import logging.config
import os
import sys
import yaml

import ROOT

import macrohelper

logger = logging.getLogger(__name__)

#______________________________________________________________________________
def draw(name, nseg=0, adcdiv=None, adcrange=None,
         tdcdiv=None, tdcrange=None, trailingdiv=None,
         totdiv=None, totrange=None, ud=True, ploop=True):
  logger.info(f'name={name}, nseg={nseg}, adcdiv={adcdiv}, '
              + f'adcrange={adcrange} '
              + f'tdcdiv={tdcdiv}, tdcrange={tdcrange}, totdiv={totdiv}, '
              + f'ud={ud}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  particle = ['', '_Pi', '_K', '_P']
  pcolor = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kRed+2]
  hitmulti = ['_OR', '_AND'] if ud else ['']
  ud = ['U', 'D'] if ud else ['']

  if adcdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(adcdiv[0], adcdiv[1])
      for i in range(nseg):
        c1.cd(i+1).SetLogy()
        for j, p in enumerate(particle):
          h1 = ROOT.gFile.Get(name + f'_ADC_seg{i}{s}{p}')
          if h1:
            if adcrange is not None:
              h1.GetXaxis().SetRangeUser(adcrange[0], adcrange[1])
            if ploop:
              h1.SetLineColor(pcolor[j])
            h1.Draw('same')
          # h2 = ROOT.gFile.Get(name + f'_AWT_seg{i}{s}{p}')
          # if h2:
          #   h2.SetLineColor(ROOT.kRed+1)
          #   h2.Draw('same')
      c1.Print(fig_path)

  if tdcdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(tdcdiv[0], tdcdiv[1])
      for i in range(nseg):
        c1.cd(i+1).SetLogy()
        for j, p in enumerate(particle):
          h1 = ROOT.gFile.Get(name + f'_TDC_seg{i}{s}{p}')
          if h1:
            if h1.GetXaxis().GetXmax() > 2000:
              h1.RebinX(4)
            if tdcrange is not None:
              h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
            if ploop:
              h1.SetLineColor(pcolor[j])
            h1.Draw('same')
      c1.Print(fig_path)

  if totdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(totdiv[0], totdiv[1])
      for i in range(nseg):
        c1.cd(i+1)
        for j, p in enumerate(particle):
          h1 = ROOT.gFile.Get(name + f'_TOT_seg{i}{s}{p}')
          if h1:
            if totrange is not None:
              h1.GetXaxis().SetRangeUser(totrange[0], totrange[1])
            if ploop:
              h1.SetLineColor(pcolor[j])
            h1.Draw('same')
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
      macrohelper.efficiency(h1)
    i = i + 1
  c1.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  macrohelper.initialize(run_info)
  draw('TriggerFlag', nseg=32, tdcdiv=(8, 4), tdcrange=(800, 1200),
       ud=False, ploop=False)
  draw('BHT', nseg=63, tdcdiv=(8, 8), totdiv=(8, 8),
       tdcrange=(1.22e6, 1.26e6), totrange=(0, 25e3))
  draw('AC', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2), tdcrange=(800, 1200),
       ud=False)
  draw('T1', nseg=1, adcdiv=(2, 2), adcrange=(0, 2000),
       tdcdiv=(2, 2), tdcrange=(1.17e6, 1.21e6))
  draw('T0', nseg=5, adcdiv=(3, 2), adcrange=(0, 2000),
       tdcdiv=(3, 2), tdcrange=(1.20e6, 1.24e6))
  draw('T0new', nseg=5, adcdiv=(3, 2), adcrange=(0, 2000),
       tdcdiv=(3, 2), tdcrange=(1.18e6, 1.22e6))
  draw('DEF', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2), tdcrange=(1.20e6, 1.25e6))
  draw('Veto', nseg=4, adcdiv=(2, 2), adcrange=(0, 2000),
       tdcdiv=(2, 2), tdcrange=(1.15e6, 1.20e6))
  draw('BTC', nseg=4, adcdiv=(2, 2), adcrange=(0, 2000),
       tdcdiv=(2, 2), tdcrange=(1.15e6, 1.20e6))
  draw('CVC', nseg=9, adcdiv=(3, 3), tdcdiv=(3, 3), tdcrange=(5e5, 7e5))
  draw('NC', nseg=6, adcdiv=(3, 2), tdcdiv=(3, 2), tdcrange=(5e5, 7e5))
  macrohelper.finalize()

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(os.path.dirname(__file__), 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  macrohelper.run(parsed.run_list, single_run)

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

mh.beamflag = [''] #

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
  hitmulti = ['_OR', '_AND'] if ud else ['']
  ud = ['U', 'D'] if ud else ['']
  if name == 'HTOF':
    ud = ['U', 'D', 'S']
  if name == 'KVC':
    ud = ['a', 'b', 'c', 'd', 'S']
  if adcdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(adcdiv[0], adcdiv[1])
      for i in range(nseg):
        c1.cd(i+1).SetLogy()
        for j, b in enumerate(mh.beamflag):
          h1 = ROOT.gFile.Get(name + f'_ADC_seg{i}{s}{b}')
          if h1:
            if adcrange is not None:
              h1.GetXaxis().SetRangeUser(adcrange[0], adcrange[1])
            if ploop:
              h1.SetLineColor(mh.beamcolor[j])
            h1.Draw('same')
          h2 = ROOT.gFile.Get(name + f'_AwT_seg{i}{s}{b}')
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
        for j, b in enumerate(mh.beamflag):
          h1 = ROOT.gFile.Get(name + f'_TDC_seg{i}{s}{b}')
          if h1:
            if h1.GetXaxis().GetXmax() > 2000:
              h1.RebinX(4)
            if tdcrange is not None:
              h1.GetXaxis().SetRangeUser(tdcrange[0], tdcrange[1])
            if ploop:
              h1.SetLineColor(mh.beamcolor[j])
            h1.Draw('same')
      c1.Print(fig_path)

  if totdiv is not None:
    for s in ud:
      c1.Clear()
      c1.Divide(totdiv[0], totdiv[1])
      for i in range(nseg):
        c1.cd(i+1)
        for j, b in enumerate(mh.beamflag):
          h1 = ROOT.gFile.Get(name + f'_TOT_seg{i}{s}{b}')
          if h1:
            if totrange is not None:
              h1.GetXaxis().SetRangeUser(totrange[0], totrange[1])
            if ploop:
              h1.SetLineColor(mh.beamcolor[j])
            h1.Draw('same')
      c1.Print(fig_path)

  c1.Clear()
  c1.Divide(2, 2)
  i = 1
  for s in hitmulti:
    c1.cd(i)
    for j, b in enumerate(mh.beamflag):
      h1 = ROOT.gFile.Get(f'{name}_HitPat{s}{b}')
      if h1:
        if ploop:
          h1.SetLineColor(mh.beamcolor[j])
        h1.Draw('same')
    i = i + 1
    c1.cd(i)
    for j, b in enumerate(mh.beamflag):
      h1 = ROOT.gFile.Get(f'{name}_Multi{s}{b}')
      if h1:
        if ploop:
          h1.SetLineColor(mh.beamcolor[j])
        h1.Draw('same')
        mh.efficiency(h1)
    i = i + 1
  c1.Print(fig_path)

#______________________________________________________________________________
def single_run(run_info):
  if os.path.basename(run_info['bin']) != 'Hodoscope':
    logger.error(f'bin must be Hodoscope: run_info={run_info}')
    return
  # comment = (f'#color[{mh.beamcolor[1]}]{{Red}}=Pion, '+
  #            f'#color[{mh.beamcolor[2]}]{{Blue}}=Kaon')
  comment = (f'#color[{ROOT.kBlack}]{{Black}}=Raw, '+
             f'#color[{ROOT.kRed+1}]{{Red}}=ADCwTDC')
  mh.initialize(run_info, __file__, comment=comment)
  draw('TriggerFlag', nseg=32, tdcdiv=(8, 4), tdcrange=(0e6, 1e6),
       ud=False, ploop=False)
  draw('BHT', nseg=63, tdcdiv=(8, 8), totdiv=(8, 8),
       tdcrange=(1.2e6, 1.6e6), totrange=(0, 25e3))
  draw('BH2', nseg=15, adcdiv=(5, 3), adcrange=(0, 2000),
       tdcdiv=(5, 3), tdcrange=(0.60e6, 0.80e6))
  draw('BAC', nseg=5, adcdiv=(3, 2), tdcdiv=(3, 2), tdcrange=(0.6e6, 0.8e6),
       ud=False)
  draw('HTOF', nseg=34, adcdiv=(6, 6), adcrange=(0, 4000),
       tdcdiv=(6, 6), tdcrange=(0.60e6, 0.80e6))
  draw('KVC', nseg=8, adcdiv=(4, 2), adcrange=(0, 4000),
       tdcdiv=(4, 2), tdcrange=(0.60e6, 0.80e6))
  draw('T1', nseg=1, adcdiv=(2, 2), adcrange=(0, 2000),
       tdcdiv=(2, 2), tdcrange=(0.60e6, 0.80e6))
  draw('CVC', nseg=8, adcdiv=(4, 2), adcrange=(0, 4000),
       tdcdiv=(4, 2), tdcrange=(0.20e6, 0.60e6))
  draw('SAC3', nseg=1, adcdiv=(2, 2), adcrange=(0, 4000),
       tdcdiv=(2, 2), tdcrange=(0.20e6, 0.60e6))
  draw('SFV', nseg=1, adcdiv=(4, 2), adcrange=(0, 4000),
       tdcdiv=(4, 2), tdcrange=(0.20e6, 0.60e6))
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

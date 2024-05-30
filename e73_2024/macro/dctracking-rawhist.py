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

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def draw_one_data(c1, name, key, nlayer):
  for i in range(nlayer):
    c1.cd(i+1) #.SetLogy()
    for j, b in enumerate(mh.beamflag):
      h1 = mh.get(f'{name}_{key}_layer{i}{b}')
      if h1:
        h1.SetLineColor(mh.beamcolor[j])
        h1.Draw('same')
        if 'TDC' in key:
          h1.GetXaxis().SetRangeUser(1200, 1500)
        if 'Multi' in key:
          mh.efficiency(h1)

#______________________________________________________________________________
@mh.update_canvas(divisions=(4, 2))
def draw_one_data2d(c1, name, key, nlayer):
  for i in range(nlayer):
    c1.cd(i+1).SetLogz()
    hname = (f'{name}_{key}_vs_HitPat_'+
             f'layer{i}{mh.beamflag_for_param}')
    h1 = mh.get(hname)
    if h1:
      h1.Draw('colz')

#______________________________________________________________________________
def draw(name, nlayer=0):
  logger.info(f'name={name}, nlayer={nlayer}')
  c1 = ROOT.gROOT.GetListOfCanvases()[0]
  fig_path = c1.GetTitle()
  for key in ['TDC', # 'Trailing',
                'TOT', 'HitPat', 'Multi']:
    for c in ['', 'C']:
      draw_one_data(name, key=c+key, nlayer=nlayer)
      if c == 'C' and (key == 'TDC' or key == 'TOT'):
        draw_one_data2d(name, key=c+key, nlayer=nlayer)

#______________________________________________________________________________
def single_run(run_info):
  comment = (f'#color[{mh.beamcolor[1]}]{{Red}}=Pion, '+
             f'#color[{mh.beamcolor[2]}]{{Blue}}=Kaon')
  mh.initialize(run_info, __file__, comment=comment)
  draw('BLC1a', nlayer=8)
  draw('BLC1b', nlayer=8)
  draw('BLC2a', nlayer=8)
  draw('BLC2b', nlayer=8)
  draw('BPC1', nlayer=8)
  draw('BPC2', nlayer=8)
  # draw('VFT', nlayer=14, tdcdiv=(5, 3))
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

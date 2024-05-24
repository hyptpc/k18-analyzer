#!/usr/bin/env python3

import logging
import multiprocessing as mp
import os
import sys

import ROOT

logger = logging.getLogger('__main__').getChild(__name__)

#______________________________________________________________________________
def fit(h1, params=None, limits=None, fitrange=(-2, 2), autozoom=False):
  f1 = ROOT.TF1('f1', '[0]*TMath::Freq((x-[1])/[2])+[3]')
  if params is not None:
    f1.SetParameters(params)
  if limits is not None:
    for i, l in enumerate(limits):
      f1.SetParLimits(i, l[0], l[1], l[2], l[3])
  for j in range(1):
    mean = f1.GetParameter(1)
    sigma = max(20, f1.GetParameter(2))
    h1.Fit('f1', '', '', mean + fitrange[0]*sigma, mean + fitrange[1]*sigma)
  if autozoom:
    h1.GetXaxis().SetRangeUser(mean-10*sigma, mean+20*sigma)
  return f1

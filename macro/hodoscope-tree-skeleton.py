#!/usr/bin/env python3

import argparse
import logging
import logging.config
import numpy as np
import os
import sys
import yaml

logger = logging.getLogger(__name__)

#______________________________________________________________________________
def run(root_path, method):
  if method == 'pyroot':
    import ROOT
    f1 = ROOT.TFile.Open(root_path)
    if f1 == None or not f1.IsOpen():
      logger.error(f'cannot find {root_path}')
      return
    tree = f1.Get('hodo')
    tree.Print()
    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      for seg in tree.t0_raw_seg:
        au = tree.t0_adc_u[seg]
        ad = tree.t0_adc_d[seg]
        tu_vec = tree.t0_tdc_u[seg]
        td_vec = tree.t0_tdc_d[seg]
        # for t in tu_vec:
        #   pass
        # for t in td_vec:
        #   pass
        logger.info(f'seg={seg}, adc_u={au}, adc_d={ad}, '
                    +f'tdc_u={tu_vec}, tdc_d={td_vec}')
    ROOT.TPython.Prompt()

  elif method == 'uproot':
    import uproot
    import matplotlib.pyplot as plt
    tree = uproot.open(root_path + ':hodo')
    ev = tree['event_number'].array()
    au = tree['t0_adc_u'].array()
    tu = tree['t0_tdc_u'].array()
    seg = 2
    # plt.hist(au[:, seg]) # all events
    # plt.show()
    # plt.plot(ev, au[:, seg]) # all events
    # plt.show()
    logger.info(tu[:,][:,])
    # plt.plot(ev, tu[:, seg])
    # plt.show()
    ''' event loop '''
    for i, seg in enumerate(tree['t0_raw_seg'].array()):
      pass
      # for j in seg:
      #   logger.info(au[i][j])
      #   logger.info(tu[i][j])


#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('root_path', help='input root path')
  parser.add_argument('method', help='method', choices=['pyroot', 'uproot'])
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(os.path.dirname(__file__), 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  run(parsed.root_path, parsed.method)

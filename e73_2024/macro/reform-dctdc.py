#!/usr/bn/env python

import argparse
import os
import logging
import sys

logger = logging.getLogger('__main__').getChild(__name__)
tmp_dir = os.path.join(os.path.dirname(__file__), 'tmp')

#______________________________________________________________________________
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('input_dctdc_path')
  parsed, unparsed = parser.parse_known_args()
  dctdc_path = parsed.input_dctdc_path #sys.argv[1]
  output_path = os.path.join(tmp_dir, os.path.basename(dctdc_path))
  with open(dctdc_path, 'r') as f:
    rows = [line for line in f]
  with open(output_path, 'w') as f:
    for line in rows:
      r = line.split()
      if len(r[0]) == 0 or r[0][0] == '#':
        if 'Cid\tPlid' in line:
          line = line.replace('Cid\tPlid', 'Layer')
        f.write(line)
      elif len(r) == 5:
        cid = int(r[0])
        plid = int(r[1])
        if cid == 101 or cid == 102:
          lid = (cid-101)*8 + plid + 1
        elif cid == 103 or cid == 104:
          lid = (cid-103)*8 + plid + 1001
        elif cid == 105 or cid == 106:
          lid = (cid-105)*8 + plid + 2001
        else:
          lid = 0
        newline = '\t'.join([str(lid), r[2], r[3], r[4]])
        f.write(newline + '\n')
    logger.info(f'generate {output_path}')

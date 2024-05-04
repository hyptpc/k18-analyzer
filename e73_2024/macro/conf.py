#!/usr/bin/env python3

import configparser
import logging
import os

logger = logging.getLogger('__main__').getChild(__name__)

#______________________________________________________________________________
def get(run_info, key=None):
  cfg = configparser.ConfigParser()
  conf_path = run_info['conf']
  logger.debug(f'read {conf_path}')
  with open(conf_path, 'r') as f:
    cfg.read_string('[DEFAULT]\n' + f.read())
  if key is None:
    return dict(cfg['DEFAULT'])
  else:
    return cfg['DEFAULT'][key.lower()]

#______________________________________________________________________________
def replace(run_info, key, value):
  conf_path = run_info['conf']
  with open(conf_path, 'r') as f:
    buf = [line for line in f]
  with open(conf_path, 'w') as f:
    for line in buf:
      sline = line.split()
      if len(sline) != 2 or sline[0][0] == '#':
        continue
      if sline[0].replace(':', '') == key:
        param_dir = os.path.dirname(sline[1])
        param_path = os.path.join(param_dir, value)
        newline = f'{key}:\t\t{param_path}'
        logger.debug(newline)
      else:
        newline = f'{sline[0]}\t\t{sline[1]}'
      f.write(newline + '\n')

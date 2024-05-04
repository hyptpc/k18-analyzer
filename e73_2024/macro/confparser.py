#!/usr/bin/env python3

import configparser
import logging

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

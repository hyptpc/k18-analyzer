#!/usr/bin/env python3

import datetime
import logging
import shutil
import os

import ROOT

import conf

logger = logging.getLogger('__main__').getChild(__name__)
tmp_dir = os.path.join(os.path.dirname(__file__), 'tmp')
os.makedirs(tmp_dir, exist_ok=True)

#______________________________________________________________________________
def output_result(run_info, result_dict, is_hrtdc=True, update=False):
  logger.debug(run_info)
  hdprm_path = conf.get(run_info, 'HDPRM')
  with open(hdprm_path, 'r') as f:
    rows = [line.split() for line in f]
  output_path = os.path.join(tmp_dir, f'HodoParam_{run_info["key"]:05d}')
  if is_hrtdc:
    p1 = -0.0009390020
  else:
    p1 = -0.8333
  comment = (f'# {datetime.datetime.now()} generated by '
             + f'{result_dict["generator"]}'
             + f' using {os.path.basename(run_info["root"])}')
  duplicate = dict()
  with open(output_path, 'w') as f:
    f.write(comment + '\n')
    for r in rows:
      if not r[0].isdigit():
        if r[0][0] != '#':
          r[0] = '#' + r[0]
        if result_dict['generator'] not in r[5]:
          f.write('\t'.join(r) + '\n')
        continue
      key = (int(r[0]), int(r[1]), int(r[2]), int(r[3]), int(r[4]))
      if key in result_dict and key not in duplicate:
        p0 = result_dict[key].GetParameter(1)
        newline = '\t'.join(r[:5] + [str(p0), str(p1)])
        logger.debug(f'{float(r[5]), float(r[6])} -> {[p0, p1]}')
        f.write(newline + '\n')
        duplicate[key] = True
      else:
        f.write('\t'.join(r) + '\n')
    logger.info(f'generate {output_path}')
  if update:
    hdprm_dir = os.path.dirname(hdprm_path)
    logger.info(f'update {os.path.join(hdprm_dir, os.path.basename(output_path))}')
    shutil.copy2(output_path, hdprm_dir)
    conf.replace(run_info, 'HDPRM', f'HodoParam_{run_info["key"]:05d}')

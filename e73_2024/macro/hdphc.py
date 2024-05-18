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
def output_result(run_info, result_dict, update=False):
  logger.debug(run_info)
  logger.debug(result_dict)
  hdphc_path = conf.get(run_info, 'HDPHC')
  with open(hdphc_path, 'r') as f:
    rows = [line.split() for line in f]
  output_path = os.path.join(tmp_dir, f'HodoPHCParam_{run_info["key"]:05d}')
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
          delimiter = ' ' if len(r[0]) == 1 else '\t'
          f.write(delimiter.join(r) + '\n')
        continue
      key = (int(r[0]), int(r[1]), int(r[2]), int(r[3]))
      if key in result_dict and key not in duplicate:
        Type = result_dict[key][0]
        nParam = result_dict[key][1]
        p0 = result_dict[key][2]
        p1 = result_dict[key][3]
        p2 = result_dict[key][4]
        newline = '\t'.join(
          r[:4] + [str(Type), str(nParam), str(p0), str(p1), str(p2)])
        if (float(r[6]), float(r[6]), float(r[7])) != (p0, p1, p2):
          logger.debug(f'update {key} {float(r[5]), float(r[6]), float(r[7])}'
                       +f'-> {(p0, p1, p2)}')
        f.write(newline + '\n')
        duplicate[key] = True
      else:
        f.write('\t'.join(r) + '\n')
    logger.info(f'generate {output_path}')
  if update:
    hdphc_dir = os.path.dirname(hdphc_path)
    logger.info(f'update {os.path.join(hdphc_dir, os.path.basename(output_path))}')
    shutil.copy2(output_path, hdphc_dir)
    conf.replace(run_info, 'HDPHC', f'HodoPHCParam_{run_info["key"]:05d}')

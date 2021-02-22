#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import copy
import logging
import os
import yaml

import classimpl
import utility

logger = logging.getLogger('__main__').getChild(__name__)

#______________________________________________________________________________
class RunlistManager(metaclass=classimpl.Singleton):
  ''' Manager class handling run list. '''

  #____________________________________________________________________________
  def __init__(self):
    self.__basename = None
    self.__work_dir = None
    self.__keys = list()
    self.__runlist = list()
    self.__is_ready = None

  #____________________________________________________________________________
  def get_bin_path(self, index):
    ''' Get value of "bin" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['bin'])
    else:
      return None

  #____________________________________________________________________________
  def get_buff_path(self, index):
    ''' Get value of "buff" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['buff'])
    else:
      return None

  #____________________________________________________________________________
  def get_conf_path(self, index):
    ''' Get value of "conf" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['conf'])
    else:
      return None

  #____________________________________________________________________________
  def get_data_path(self, index):
    ''' Get value of "data" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['data'])
    else:
      return None

  #____________________________________________________________________________
  def get_div_unit(self, index):
    ''' Get value of "unit" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return int(self.__runlist[index]['unit'])
    else:
      return None

  #____________________________________________________________________________
  def get_key(self, index):
    ''' Get value of "key" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['key'])
    else:
      return None

  #____________________________________________________________________________
  def get_keys(self):
    ''' Get list of keys '''
    return list(self.__keys)

  #____________________________________________________________________________
  def get_nevents(self, index):
    ''' Get value of "nevents" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return int(self.__runlist[index]['nevents'])
    else:
      return None

  #____________________________________________________________________________
  def get_nevents_recorder(self, path, runno=None):
    ''' Get number of events of raw data in recorder.log. '''
    if self.__is_ready is False:
      return None
    nevents = None
    if os.path.isdir(path):
      reclog_path = os.path.join(path, 'recorder.log')
      if os.path.isfile(reclog_path):
        cand = list()
        with open(reclog_path, 'r') as freclog:
          for line in freclog:
            words = line.split()
            if len(words) > 2 and runno == int(words[1]):
              cand.append(words[15]) if len(words) > 15 else -1
        nevents = int(cand[0]) if len(cand) == 1 else None
    return nevents

  #____________________________________________________________________________
  def get_nproc(self, index):
    ''' Get value of "nproc" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return int(self.__runlist[index]['nproc'])
    else:
      return None

  #____________________________________________________________________________
  def get_nruns(self):
    ''' Get length of keys '''
    return len(self.__keys)

  #____________________________________________________________________________
  def get_queue(self, index):
    ''' Get value of "queue" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['queue'])
    else:
      return None

  #____________________________________________________________________________
  def get_root_path(self, index):
    ''' Get value of "root" in run list. '''
    index = self.__get_index(index)
    if index is not None:
      return str(self.__runlist[index]['root'])
    else:
      return None

  #____________________________________________________________________________
  def get_run_info(self, index):
    ''' Get run info dictionary of run list. '''
    index = self.__get_index(index)
    if index is not None:
      return dict(self.__runlist[index])
    else:
      return None

  #____________________________________________________________________________
  def get_tag(self):
    ''' Get tag (filehead of run list yaml). '''
    return self.__basename

  #____________________________________________________________________________
  def get_work_dir(self, path):
    ''' Get directory of "WORKDIR" in run list yaml. '''
    if self.__is_ready is False:
      return None
    data = dict()
    with open(path, 'r') as f:
      data = yaml.safe_load(f.read())
    work_dir = os.path.expanduser(data['WORKDIR'])
    if os.path.isdir(work_dir):
      return work_dir
    else:
      logger.error(f'Cannot find {work_dir}')
      exit(1)

  #____________________________________________________________________________
  def set_run_list(self, path):
    ''' Set run list yaml. '''
    self.__basename = os.path.splitext(os.path.basename(path))[0]
    self.__work_dir = self.get_work_dir(path)
    self.__keys = list()
    self.__runlist = list()
    self.__make_run_list(path)

  #____________________________________________________________________________
  def __decode_run_list(self, path):
    if self.__is_ready is False:
      return None
    if self.__work_dir is None:
      self.__work_dir = self.get_work_dir(path)
    data = dict()
    with open(path, 'r') as f:
      data = yaml.safe_load(f.read())
    defset = data['DEFAULT']
    runlist = list()
    for key, parsets in data['RUN'].items():
      if parsets is None:
        runlist.append([key, defset])
      else:
        runset = copy.deepcopy(defset)
        runset.update(parsets)
        runlist.append([key, runset])
    return runlist

  #____________________________________________________________________________
  def __get_index(self, index):
    ''' get index from keys. '''
    # if index in self.__keys:
    #   return self.__keys.index(index)
    # elif
    if isinstance(index, int) and len(self.__runlist) > index:
      return index
    else:
      return None

  #____________________________________________________________________________
  def __make_data_path(self, path, runno=None):
    ''' Make raw data file path. '''
    if self.__is_ready is False:
      return None
    if not os.path.exists(path):
      logger.error(f'Cannot find file: {path}')
    else:
      if os.path.isfile(path):
        return os.path.realpath(path)
      elif (os.path.isdir(path)
            and runno is not None
            and isinstance(runno, int)):
        base = f'{path}/run{runno:05d}'
        data = f'{base}.dat.gz'
        if not os.path.isfile(data):
          data = f'{base}.dat'
        if os.path.isfile(data):
          return os.path.realpath(data)
        else:
          logger.error(f'Cannot find file: {data}')
    logger.error('Cannot decide deta file path')
    exit(1)

  #____________________________________________________________________________
  def __make_dstin_path(self, base, key_array):
    ''' Make path array of input files for dst analysis. '''
    if self.__is_ready is False:
      return None
    if len(key_array) == 0:
      logger.error('Cannot decide input files for dst')
      exit(1)
    dstin_path = list()
    for key in key_array:
      split = os.path.splitext(base)
      dstin_path.append(split[0] + key + split[1])
    return dstin_path

  #____________________________________________________________________________
  def __make_root_path(self, path, base=None):
    ''' Make output root file path. '''
    if self.__is_ready is False:
      return None
    root_path = None
    if not os.path.exists(path):
      dir_path = os.path.dirname(path)
      if os.path.isdir(dir_path):
        return os.path.realpath(path)
    elif os.path.isfile(path):
      return os.path.realpath(path)
    elif os.path.isdir(path) and base is not None:
      return os.path.realpath(os.path.join(path, base+'.root'))
    logger.error('Cannot decide root file path')
    exit(1)

  #____________________________________________________________________________
  def __make_run_list(self, path):
    ''' Make run list. '''
    if self.__is_ready is False:
      return
    cdir = os.getcwd()
    raw_runlist = self.__decode_run_list(path)
    os.chdir(self.__work_dir)
    for item in raw_runlist:
      run = dict()
      run['key'] = item[0]
      pbin = None
      if os.path.isfile(item[1]['bin']):
        pbin = item[1]['bin']
      else:
        logger.error(f"Cannot find file: {item[1]['bin']}")
        exit(1)
      run['bin'] = pbin
      runno = None
      if os.path.isfile(item[1]['data']):
        tmp = os.path.splitext(os.path.basename(item[1]['data']))[0]
        runno = int(tmp[3:8]) if tmp[3:8].isdigit() else None
      else:
        runno = item[0] if isinstance(item[0], int) else None
      pdata = self.__make_data_path(item[1]['data'], runno)
      run['data'] = pdata
      run['nevents'] = self.get_nevents_recorder(
        os.path.dirname(os.path.abspath(pdata)), runno)
      pconf = None
      if os.path.exists(item[1]['conf']):
        if os.path.isfile(item[1]['conf']):
          pconf = item[1]['conf']
        elif os.path.isdir(item[1]['conf']) and runno is not None:
          pconf = os.path.join(item[1]['conf'], f'analyzer_{runno:05d}.conf')
      if pconf is None:
        logger.error('Cannot decide conf file path')
        exit(1)
      run['conf'] = pconf
      base = (item[0] + os.path.basename(pbin) if runno is None
              else f'run{runno:05d}_{os.path.basename(pbin)}')
      run['root'] = self.__make_root_path(item[1]['root'], base)
      if 'Dst' in run['bin']:
        if 'dstin' in item[1]:
          base = run['root'].replace(os.path.basename(run['bin']), '')
          run['dstin'] = self.__make_dstin_path(base, item[1]['dstin'])
        else:
          logger.error(f'{run["bin"]} needs input files set as "dstin".')
          exit(1)
      run['unit'] = (item[1]['unit'] if isinstance(item[1]['unit'], int)
                     else 0)
      run['queue'] = (item[1]['queue'] if isinstance(item[1]['queue'], str)
                      else 's')
      run['nproc'] = (item[1]['nproc'] if isinstance(item[1]['nproc'], int)
                      else 1)
      if item[1]['buff'] is None:
        pbuff = None
      elif os.path.isdir(item[1]['buff']):
        pbuff = item[1]['buff']
      else:
        logger.error('Cannot decide buffer file path')
        exit(1)
      run['buff'] = pbuff
      self.__runlist.append(run)
      self.__keys.append(run['key'])
    os.chdir(self.__work_dir)

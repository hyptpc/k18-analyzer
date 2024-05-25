#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import classimpl
import logging
import shlex
import subprocess
import threading
import time

BSUB_RESPONCE = 'Job <JobID> is submitted to queue <queue>'
logger = logging.getLogger('__main__').getChild(__name__)

#______________________________________________________________________________
class BJobManager(metaclass=classimpl.Singleton):
  ''' BJob class throwing subprocess of "bjobs" and "bkill". '''

  #____________________________________________________________________________
  def __init__(self):
    self.__interval = 10 # [s]
    self.__updater_thread = threading.Thread(target=self.__updater)
    self.__updater_thread.name = 'updater'
    self.__updater_thread.daemon = True
    self.__updater_status = 'IDLE'
    self.__buf = None
    self.__status_list = dict()

  #______________________________________________________________________________
  def __del__(self):
    while self.__updater_thread.is_alive():
      self.__updater_thread.join()
    logger.debug('bye')

  #____________________________________________________________________________
  def __updater(self):
    logger.debug('updater thread starts')
    ptime = time.time()
    while self.__updater_status == 'RUNNING':
      logger.debug(f'updater is running')
      self.update_job_status()
      sleep = self.__interval - (time.time() - ptime)
      if sleep > 0:
        time.sleep(sleep)
      ptime = time.time()
    logger.debug('updater thread ends')

  #______________________________________________________________________________
  def stop(self):
    ''' Stop updater thread. '''
    self.updater_status = 'END'

  #______________________________________________________________________________
  def get_job_status(self, job_id, start_time):
    ''' Get job status. '''
    while self.__updater_status == 'UPDATING':
      time.sleep(0.2)
    if self.__updater_status == 'RUNNING':
      s = self.__status_list.get(job_id, None)
      if s is not None:
        return s
      if time.time() - start_time > 3000:
        logger.debug(f'{job_id} is missing')
        return 'DONE'
    return 'INIT'

  #______________________________________________________________________________
  def update_job_status(self):
    ''' Update job status. '''
    self.__updater_status = 'UPDATING'
    logger.debug('updating bjob status')
    cmd = 'bjobs -a'
    proc = None
    try:
      proc = subprocess.run(shlex.split(cmd),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=True)
    except subprocess.CalledProcessError as e:
      logger.error(e)
      self.__status = 'FAILED'
      return self.__status
    self.__buf = proc.stdout.decode()
    self.__status_list = dict()
    for line in self.__buf.splitlines():
      columns = line.split()
      if len(columns) < 4:
        continue
      if not columns[0].isdigit():
        continue
      job_id = int(columns[0])
      self.__status_list[job_id] = columns[2]
    logger.debug('updated bjob status')
    self.__updater_status = 'RUNNING'

  #____________________________________________________________________________
  def run(self):
    self.__updater_status = 'RUNNING'
    self.__updater_thread.start()

  #____________________________________________________________________________
  def start(self):
    self.run()

  #____________________________________________________________________________
  @staticmethod
  def read_job_id(buff):
    ''' Read job id. '''
    job_id = None
    flag = True
    words = buff.split()
    for i, item in enumerate(BSUB_RESPONCE.split()):
        if i == 1 or i == 6:
          continue
        if item != words[i]:
          flag = False
    if flag is True:
      job_id = int(words[1][1:-1])
    return job_id

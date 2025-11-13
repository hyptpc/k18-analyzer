#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import argparse
import fcntl
import json
import logging
import logging.config
import os
import sys
import time
import yaml
from datetime import datetime, timedelta

top_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(top_dir, 'module'))

import runmanager
import singlerun
import utility
from utility import pycolor as cl
import notify

DISPLAY_PERIOD = 5  # second
logger = logging.getLogger(__name__)

#______________________________________________________________________________
def display(filename, webhook_url, run_name, last_statuses, use_notify, initial_run, first_loop_completed):
  buff = str()
  with open(filename, 'r') as f:
    try:
      fcntl.flock(f.fileno(), fcntl.LOCK_SH)
    except IOError as e:
      logger.error(e)
      return False, {}, {}
    else:
      buff = f.read()
    finally:
      fcntl.flock(f.fileno(), fcntl.LOCK_UN)
  try:
    info = json.loads(buff)
  except ValueError as e:
    logger.error(e)
    return False, {}, {}
  os.system('clear')
  buff = (cl.reverce + cl.bold
          + 'KEY'.ljust(8) + '  '
          + 'STAT'.ljust(18) + '  '
          + 'BIN'.ljust(16) + '  '
          + 'CONF'.ljust(16) + '  '
          + 'DATA(#EVENT)'.ljust(24) + '  '
          + 'ROOT'.ljust(16) + '  '
          + 'TIME'.ljust(8)
          + cl.end)
  print(buff, flush=True)
  
  n_unfinished = 0
  status_counts = {}
  current_statuses = {}
  
  for key, item in sorted(info.items(),
                          key=lambda x:int(x[0]) if isinstance(x[0], int) else x[0]):
    status_str = singlerun.SingleRun.decode_status(item)
    current_statuses[key] = status_str

    if use_notify:
      last_status = last_statuses.get(key)
      if status_str != last_status:
        if not (initial_run and not first_loop_completed):
            current_main_status = status_str.split(':')[0].split('(')[0]
            last_main_status = ""
            if last_status is not None:
                last_main_status = last_status.split(':')[0].split('(')[0]

            if not (last_main_status == "running" and current_main_status == "running"):
                if last_status is None and 'staged' in status_str:
                  notify.send_to_discord(webhook_url,
                    f"[ JOB SUBMITTED ] batch '{run_name}': `{key}`")
                elif 'error' in status_str or 'terminated' in status_str:
                  notify.send_to_discord(webhook_url,
                    f"[ JOB FAILED ] batch '{run_name}': `{key}` is now `{status_str}`")
                elif last_status is not None:
                  notify.send_to_discord(webhook_url,
                    f"[ STATUS CHANGE ] batch '{run_name}': `{key}`  {last_status} → {status_str}")
                    
    simple_status = status_str
    if simple_status.startswith('running'):
      simple_status = 'running'
    elif simple_status.startswith('merging'):
      simple_status = 'merging'
    elif simple_status.startswith('done'):
      simple_status = 'done'
    status_counts[simple_status] = status_counts.get(simple_status, 0) + 1
          
    if 'done' not in status_str and 'error' not in status_str:
      n_unfinished += 1
    
    ptime = singlerun.SingleRun.decode_time(item)
    infile = None
    if 'data' in item and item['data'] is not None:
      infile = item['data']
    elif 'dstin' in item and len(item['dstin']) == 1:
      infile = item['dstin'][0]
    nev = None
    if 'nev' in item:
      nev = item['nev']
    inbuf = f'{os.path.basename(str(infile))} ({str(nev)})'
    buff = (cl.bold + key[:8].ljust(8) + cl.end + '  '
            + f'{cl.reverce}{cl.bold}{cl.red}{status_str}{cl.end}'.ljust(16 + 20)
            + '  '
            + os.path.basename(item['bin'])[-16:].ljust(16) + '  '
            + os.path.basename(item['conf'])[-16:].ljust(16) + '  '
            + inbuf[-24:].ljust(24)
            + '  '
            + os.path.basename(item['root'])[-16:].ljust(16) + '  '
            + ptime.rjust(8))
    print(buff,flush=True)

  total_jobs = len(info)
  summary_message = f"Progress: {total_jobs - n_unfinished} / {total_jobs} complete"
  details = [f"{stat}: {count}" for stat, count in status_counts.items()]
  summary_message += f" ({', '.join(details)})"
    
  print('\n' + cl.reverce + cl.bold + summary_message + cl.end, flush=True)
  print(cl.reverce + cl.bold + "Press 'Ctrl-C' to exit" + cl.end, flush=True)
  
  return n_unfinished, summary_message, current_statuses

#______________________________________________________________________________
def main(path, webhook_url, interval_hours, use_notify, auto_exit, initial_run):
  
  if use_notify:
      lock_file_path = path + ".lock"

      if os.path.exists(lock_file_path):
          try:
              with open(lock_file_path, 'r') as f:
                  pid_in_lock = int(f.read().strip())
              os.kill(pid_in_lock, 0)
              logger.warning(f"Monitor for {os.path.basename(path)} with notifications is already running with PID {pid_in_lock}. Exiting.")
              return
          except (IOError, ValueError, OSError):
              logger.info("Found a stale lock file. Taking over.")
              
      try:
        with open(lock_file_path, 'w') as f:
            f.write(str(os.getpid()))

        print('monitor started')
        
        ptime = time.time()
        last_notification_time = datetime.now()
        if parsed.interval_sec is not None:
          notification_interval = timedelta(seconds=parsed.interval_sec)
        else:
          notification_interval = timedelta(hours=interval_hours)
        last_statuses = {}
        last_summary = ""
        first_loop_completed = False
        run_name = os.path.splitext(os.path.basename(path))[0]

        if use_notify and not initial_run:
          notify.send_to_discord(webhook_url,f"[ MONITORING STARTED ] Now watching analysis batch: '{run_name}'")

        try:
              while True:
                n_unfinished, summary, current_statuses = display(
                  path, webhook_url, run_name, last_statuses, use_notify, initial_run, first_loop_completed
                )
                if use_notify and n_unfinished > 0 and (datetime.now() - last_notification_time > notification_interval):
                  is_stalled = "running" in last_summary and summary == last_summary
                  if is_stalled:
                    logger.info('Job is stalled. Sending status report.')
                    # notify.send_to_discord(webhook_url, f"[ STATUS REPORT ] For '{run_name}':\n> {summary}")
                  last_notification_time = datetime.now()
                last_statuses = current_statuses
                last_summary = summary
                if n_unfinished == 0 and len(last_statuses) > 0:
                  if auto_exit:
                    logger.info('All jobs have been completed. Exiting automatically.')
                    if use_notify and not (initial_run and not first_loop_completed):
                      notify.send_to_discord(webhook_url, f"**[ BATCH COMPLETE ] All analysis jobs for '{run_name}' have completed!**")
                    break
                  else:
                    print(cl.reverce + cl.bold + "All jobs are complete. Press 'Ctrl-C' to exit." + cl.end, flush=True)
                first_loop_completed = True
                dtime = DISPLAY_PERIOD - (time.time() - ptime)
                if dtime > 0:
                  time.sleep(dtime)
                ptime = time.time()
        except KeyboardInterrupt as e:
          logger.info(e)
          if use_notify:
            notify.send_to_discord(webhook_url, f"[ MONITORING STOPPED ] Monitoring for '{run_name}' was stopped manually.")
        except FileNotFoundError as e:
          logger.error(e)
          if use_notify:
            notify.send_to_discord(webhook_url, f"[ ERROR ] Status file not found for '{run_name}': ({path})")
        finally:
          print(f"\nStopping monitor for: {run_name}")

      finally:
        if os.path.exists(lock_file_path):
            os.remove(lock_file_path)

  else:
      print('monitor started')
      
      ptime = time.time()
      last_notification_time = datetime.now()
      if parsed.interval_sec is not None:
        notification_interval = timedelta(seconds=parsed.interval_sec)
      else:
        notification_interval = timedelta(hours=interval_hours)
      last_statuses = {}
      last_summary = ""
      first_loop_completed = False
      run_name = os.path.splitext(os.path.basename(path))[0]

      if use_notify and not initial_run:
        notify.send_to_discord(webhook_url,f"[ MONITORING STARTED ] Now watching analysis batch: '{run_name}'")

      try:
            while True:
              n_unfinished, summary, current_statuses = display(
                path, webhook_url, run_name, last_statuses, use_notify, initial_run, first_loop_completed
              )
              if use_notify and n_unfinished > 0 and (datetime.now() - last_notification_time > notification_interval):
                is_stalled = "running" in last_summary and summary == last_summary
                if is_stalled:
                  logger.info('Job is stalled. Sending status report.')
                  # notify.send_to_discord(webhook_url, f"[ STATUS REPORT ] For '{run_name}':\n> {summary}")
                last_notification_time = datetime.now()
              last_statuses = current_statuses
              last_summary = summary
              if n_unfinished == 0 and len(last_statuses) > 0:
                if auto_exit:
                  logger.info('All jobs have been completed. Exiting automatically.')
                  if use_notify and not (initial_run and not first_loop_completed):
                    notify.send_to_discord(webhook_url, f"**[ BATCH COMPLETE ] All analysis jobs for '{run_name}' have completed!**")
                  break
                else:
                  print(cl.reverce + cl.bold + "All jobs are complete. Press 'Ctrl-C' to exit." + cl.end, flush=True)
              first_loop_completed = True
              dtime = DISPLAY_PERIOD - (time.time() - ptime)
              if dtime > 0:
                time.sleep(dtime)
              ptime = time.time()
      except KeyboardInterrupt as e:
        logger.info(e)
        if use_notify:
          notify.send_to_discord(webhook_url, f"[ MONITORING STOPPED ] Monitoring for '{run_name}' was stopped manually.")
      except FileNotFoundError as e:
        logger.error(e)
        if use_notify:
          notify.send_to_discord(webhook_url, f"[ ERROR ] Status file not found for '{run_name}': ({path})")
      finally:
        print(f"\nStopping monitor for: {run_name}")

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('stat_json', help='stat json file path')
  parser.add_argument('--webhook', default=None, help='Discord webhook URL. Required for notifications.')
  parser.add_argument('--interval', type=int, default=12, help='Notification interval in hours for progress reports.')
  parser.add_argument('--interval-sec', type=float, default=None, help='Notification interval in seconds (overrides --interval if set).')
  parser.add_argument('--notify', action='store_true', help='Enable Discord notifications for this monitor instance.')
  parser.add_argument('--auto-exit', action='store_true', help='Exit automatically when all jobs are complete.')
  parser.add_argument('--initial-run', action='store_true',
                      help='Indicates an initial run for a pre-existing file to suppress certain notifications.')
  
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(top_dir, 'logging_config.yml')
  if os.path.exists(log_conf):
    with open(log_conf, 'r') as f:
      logging.config.dictConfig(yaml.safe_load(f))
  else:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

  # Enable notifications only if both --notify flag and --webhook URL are provided.
  notifications_enabled = parsed.notify and parsed.webhook is not None
  if parsed.notify and parsed.webhook is None:
    logger.warning("--notify flag was provided, but --webhook URL is missing. No notifications will be sent.")

  main(parsed.stat_json, parsed.webhook, parsed.interval, notifications_enabled, parsed.auto_exit, parsed.initial_run)

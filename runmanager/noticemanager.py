#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__date__ = '05 Oct. 2025'

import os
import sys
import argparse
import json
import time
import subprocess
import signal
import atexit

# --- Configuration ---
APP_DIR = os.path.dirname(os.path.abspath(__file__))
STAT_DIR = os.path.join(APP_DIR, 'stat')
MONITOR_SCRIPT_PATH = os.path.join(APP_DIR, 'monitor.py')
PYTHON_EXECUTABLE = sys.executable
PID_FILE = os.path.join(APP_DIR, '.noticemanager.pid')
KEYWORD_FILE_PATH = os.path.join(APP_DIR, 'keywords.txt')

sys.path.append(os.path.join(APP_DIR, 'module'))
import singlerun

def is_process_running(pid):
    try: os.kill(pid, 0)
    except OSError: return False
    else: return True

def is_batch_already_completed(file_path, log_file_handle):
    try:
        with open(file_path, 'r') as f: data = json.load(f)
        if not data: return True
        for job_info in data.values():
            status = singlerun.SingleRun.decode_status(job_info)
            if 'done' not in status and 'terminated' not in status:
                return False
        return True
    except (json.JSONDecodeError, IOError, KeyError) as e:
        log_file_handle.write(f"    [WARNING] Could not parse {os.path.basename(file_path)}. Will attempt to monitor. Error: {e}\n")
        return False

def cleanup():
    if os.path.exists(PID_FILE): os.remove(PID_FILE)
    print("Daemon stopped and PID file cleaned up.")

def start_manager(args):
    if os.path.exists(PID_FILE):
        with open(PID_FILE, 'r') as f: pid = int(f.read())
        if is_process_running(pid):
            print(f"Notice Manager is already running with PID: {pid}")
            return

    keywords_to_monitor, monitoring_mode = [], ""
    if args.use_keyword_file:
        if os.path.exists(KEYWORD_FILE_PATH):
            monitoring_mode = f"file ({os.path.basename(KEYWORD_FILE_PATH)})"
            with open(KEYWORD_FILE_PATH, 'r') as f:
                for line in f:
                    stripped_line = line.strip()
                    if stripped_line and not stripped_line.startswith('#'):
                        keywords_to_monitor.append(stripped_line)
        else: monitoring_mode = "all (keyword file not found)"
    else:
        monitoring_mode = f"keyword ('{args.keyword}')"
        if args.keyword.lower() != 'all': keywords_to_monitor.append(args.keyword)
    
    if os.fork() > 0: sys.exit()
    os.setsid()
    if os.fork() > 0: sys.exit()

    with open(PID_FILE, 'w') as f: f.write(str(os.getpid()))
    atexit.register(cleanup)
    
    debug_log_path = os.path.join(APP_DIR, 'noticemanager_debug.log')
    with open(debug_log_path, 'a+') as log_file:
        log_file.seek(0); log_file.truncate()
        log_file.write(f"--- Notice Manager started at {time.asctime()} (PID: {os.getpid()}) ---\n")
        log_file.write(f"    Monitoring mode: {monitoring_mode}\n")
        running_monitors = {}
        initial_files = set()
        if os.path.isdir(STAT_DIR):
            initial_files = {filename for filename in os.listdir(STAT_DIR) if filename.endswith('.json')}
        while True:
            try:
                completed_monitors = []
                for filename, proc in running_monitors.items():
                    if proc.poll() is not None:
                        completed_monitors.append(filename)
                for filename in completed_monitors:
                    del running_monitors[filename]

                if os.path.isdir(STAT_DIR):
                    for filename in os.listdir(STAT_DIR):
                        if not filename.endswith('.json'): continue
                        if keywords_to_monitor:
                            if not any(kw in filename for kw in keywords_to_monitor): continue
                        if filename in running_monitors: continue
                        stat_file_path = os.path.join(STAT_DIR, filename)
                        if is_batch_already_completed(stat_file_path, log_file): continue
                        cmd = [
                            PYTHON_EXECUTABLE, MONITOR_SCRIPT_PATH, stat_file_path,
                            '--webhook', args.webhook, '--auto-exit', '--notify'
                        ]
                        if filename in initial_files: cmd.append('--initial-run')
                        log_file.write(f"\n--- Spawning monitor for {filename} at {time.asctime()} ---\n")
                        proc = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT)
                        running_monitors[filename] = proc
            except Exception as e:
                log_file.write(f"\n!!! UNEXPECTED ERROR in main loop at {time.asctime()} !!!\n")
                log_file.write(f"    Error: {e}\n")
            time.sleep(10)

def stop_all(args):
    if not os.path.exists(PID_FILE):
        print("Notice Manager is not running.")
        os.system(f"pkill -f '{MONITOR_SCRIPT_PATH}'")
        return
    with open(PID_FILE, 'r') as f: pid = int(f.read())
    if is_process_running(pid):
        print(f"Stopping Notice Manager daemon (PID: {pid})...")
        try: os.kill(pid, signal.SIGTERM)
        except OSError: pass
    else:
        print("Notice Manager is not running (stale PID file found).")
    print("--- Stopping all monitor processes... ---")
    os.system(f"pkill -f '{MONITOR_SCRIPT_PATH}'")
    cleanup()

def show_status(args):
    if os.path.exists(PID_FILE):
         with open(PID_FILE, 'r') as f: pid = int(f.read())
         if is_process_running(pid): print(f"Notice Manager is [ Running ] (PID: {pid})")
         else: print("Notice Manager is [ Stopped ] (stale PID file found)")
    else:
        print("Notice Manager is [ Stopped ]")

    print("--- Child Monitor Process Status ---")
    monitors_found = False
    try:
        result = subprocess.run(['pgrep', '-af', MONITOR_SCRIPT_PATH], capture_output=True, text=True, check=True)
        
        for line in result.stdout.strip().split('\n'):
            if not line: continue
            parts = line.split()
            pid_str = parts[0]
            
            json_file_path = ""
            for part in parts:
                if part.endswith('.json') and STAT_DIR in part:
                    json_file_path = part
                    break
            
            if json_file_path:
                pid = int(pid_str)
                status = "[ Running ]"
                print(f"  {status} PID: {pid:<7} | File: {os.path.basename(json_file_path)}")
                monitors_found = True

    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    if not monitors_found:
        print("-> No running monitor processes found.")
# -----------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A tool to manage multiple analysis job monitors.')
    subparsers = parser.add_subparsers(dest='command', required=True)
    parser_start = subparsers.add_parser('start', help='Start the notice manager daemon.')
    parser_start.add_argument('--webhook', required=True, help='Discord webhook URL')
    parser_start.add_argument('--keyword', default='all', help='Keyword to filter files. Default: "all".')
    parser_start.add_argument('--use-keyword-file', action='store_true', help='Use the keywords.txt file for filtering.')
    parser_start.set_defaults(func=start_manager)
    parser_stop = subparsers.add_parser('stop', help='Stop the notice manager daemon and all monitors.')
    parser_stop.set_defaults(func=stop_all)
    parser_status = subparsers.add_parser('status', help='Show the status of the manager and monitors.')
    parser_status.set_defaults(func=show_status)
    args = parser.parse_args()
    args.func(args)

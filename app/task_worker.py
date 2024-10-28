#!/usr/bin/env python

import traceback
import subprocess
from time import time

from __main__ import (
  get_lock,
  task_queue,
  master_lock,
  task_status)

def task_worker():
    while True:
        task_id, task_spec, (task, data) = task_queue.get()
        with master_lock:
            lock = get_lock("task_status")
        with lock:
            task_status["specs"].append(task_spec)
            task_status[task_id] = {
              "progress": "",
              "code": 202,
              "start": time()}

        try:
            task(**data)
        except subprocess.CalledProcessError:
            with lock:
                task_status[task_id]["code"] = 500
        except Exception as e:
            with lock:
                task_status[task_id]["progress"] += "\n" + traceback.format_exc()
                task_status[task_id]["code"] = 500
        else:
            with lock:
                task_status[task_id]["code"] = 200
        finally:
            task_queue.task_done()
            task_status["specs"].remove(task_spec)

def progress(task_id, msg, bold=True, newline=True):
    if bold:
        msg = f"<b>{msg}</b>"
    if newline:
        msg += "\n"
    with master_lock:
        lock = get_lock("task_status")
    with lock:
        task_status[task_id]["progress"] += msg
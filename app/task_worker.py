#!/usr/bin/env python

import traceback
import subprocess

from __main__ import (
  lock,
  task_queue,
  task_status)

def task_worker():
    while True:
        task_id, task_spec, (task, data) = task_queue.get()
        with lock:
            task_status["specs"].append(task_spec)
            task_status[task_id] = {"progress": "", "code": 202}

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

def progress(task_id, msg, bold=True):
    if bold:
        msg = f"<b>{msg}</b>"
    msg += "\n"
    with lock:
        task_status[task_id]["progress"] += msg
#!/usr/bin/env python

import re
import os
import pty
import tty
#import traceback
import subprocess
from time import time

from flask_app import (
    get_lock,
    running_pr,
    task_queue,
    master_lock,
    task_status,
    cancelled_tasks,
)

from app.dirs import data_dir

###############################################################################
### Run a task
###############################################################################


def task_worker():
    while True:
        with master_lock:
            queue_lock = get_lock("__task_queue__")
            status_lock = get_lock("__task_status__")
        with queue_lock:
            task_id, task_spec, task_params, cmd = task_queue.get()
            genome = task_spec[1]

        if task_id in cancelled_tasks:
            with queue_lock:
                task_queue.task_done()
            with status_lock:
                if task_spec in task_status["specs"]:
                    task_status["specs"].remove(task_spec)
            continue

        with status_lock:
            task_status["specs"].append(task_spec)
            task_status[task_id] = {"progress": "", "code": 202, "start": time()}

        params_str = "\n".join([" ".join(param) for param in task_params]) + "\n\n"
        task_status[task_id]["progress"] += "<b>Task parameters</b>\n"
        task_status[task_id]["progress"] += params_str

        #        if task_spec[0] == "setup":
        #            log_root = f"{genome}.setup"
        #        elif task_spec[0] == "filt_qc":
        #            log_root = f"{genome}.{task_spec[2]}.preparation"
        #        else:
        #            log_root = task_id
        log_fn = f"{task_id}.log"

        log_file = open(os.path.join(data_dir, genome, "logs", log_fn), "a")
        _ = log_file.write(params_str)

        show_stdout = True
        master_fd, slave_fd = pty.openpty()
        tty.setraw(slave_fd)

        pr = subprocess.Popen(
            cmd,
            stdout=slave_fd,
            stderr=subprocess.STDOUT,
            preexec_fn=os.setsid,
            close_fds=True,
        )
        os.close(slave_fd)

        with master_lock:
            running_pr_lock = get_lock("__running_pr__")
        with running_pr_lock:
            running_pr[task_id] = (pr, task_spec)

        if show_stdout:
            buffer = b""
            with os.fdopen(master_fd, "rb") as master:
                while True:
                    try:
                        char = master.read(1)
                    except OSError as e:
                        if e.errno == 5:
                            break
                        else:
                            raise
                    else:
                        if not char:
                            break

                        if char in (b"\n", "\n"):
                            buffer += char
                            buffer_str = buffer.decode("utf-8", errors="replace")
                            with status_lock:
                                task_status[task_id]["progress"] += buffer_str
                            _ = log_file.write(
                                buffer_str.replace("<b>", "").replace("</b>", "")
                            )

                            buffer = b""

                        else:
                            if char in (b"\r", b"\x08"):
                                if char == b"\r":
                                    next_char = master.read(1)
                                    if next_char and next_char == b"\n":
                                        char = next_char
                                    else:
                                        char += next_char

                                buffer_str = buffer.decode("utf-8", errors="replace")
                                with status_lock:
                                    task_status[task_id]["progress"] += buffer_str
                                _ = log_file.write(
                                    buffer_str.replace("<b>", "").replace("</b>", "")
                                )

                                buffer = b""
                            buffer += char
            if buffer:
                buffer_str = buffer.decode("utf-8", errors="replace")
                # Remove ANSI escape codes
                buffer_str = re.sub(
                    "\x1b[@-_][0-?]*[ -/]*[@-~]", "", buffer_str
                ).strip()
                with status_lock:
                    task_status[task_id]["progress"] += buffer_str
                _ = log_file.write(buffer_str.replace("<b>", "").replace("</b>", ""))

        pr.wait()

        end = time()
        task_status[task_id]["end"] = end
        secs = end - task_status[task_id]["start"]

        hours = secs // 3600
        minutes = (secs % 3600) // 60
        seconds = round(secs % 60)

        duration = (
            (f"{hours} h " if hours else "")
            + (f"{minutes} min " if minutes else "")
            + f"{seconds} sec"
        )

        _ = log_file.write("\nTotal time: " + duration)
        log_file.close()

        with running_pr_lock:
            _ = running_pr.pop(task_id, None)

        # if show_stdout:
        #    with status_lock:
        #        task_status[task_id]["progress"] += "\n"

        with status_lock:
            if task_status[task_id]["code"] != 410:
                if pr.returncode == 0:
                    task_status[task_id]["code"] = 200
                else:
                    task_status[task_id]["code"] = 500

        # task_status[task_id]["progress"] += "\n" + traceback.format_exc()
        # task_status[task_id]["code"] = 500

        with queue_lock:
            task_queue.task_done()
        with status_lock:
            if task_spec in task_status["specs"]:
                task_status["specs"].remove(task_spec)

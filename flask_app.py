#!/usr/bin/env python

import queue
from flask import Flask
from threading import Thread, Lock
from flask_compress import Compress

app = Flask(__name__)
app.secret_key = "275b250b2816223102a047d297676ce0"
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024

Compress(app)

lock_dict = {}
master_lock = Lock()


def get_lock(obj):
    return lock_dict.setdefault(obj, Lock())


running_pr = dict()
task_queue = queue.Queue()
task_status = {"specs": []}
cancelled_tasks = []

from app.task_worker import task_worker

worker_thread = Thread(target=task_worker, daemon=True)
worker_thread.start()

from app.routes import (
    assoc,
    error,
    index,
    fa_res,
    fa_opts,
    assoc_qc,
    assoc_res,
    assoc_filt,
    assoc_opts,
)

if __name__ == "__main__":
    app.run(debug=True, port=5101)

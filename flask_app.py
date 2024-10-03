#!/usr/bin/env python

import os
import queue
from flask import Flask
from threading import Thread, Lock
from flask_compress import Compress

from app.dirs import (
  tmp_dir,
  data_dir,
  root_dir)

app = Flask(__name__)
app.secret_key = "275b250b2816223102a047d297676ce0"
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024
app.config["UPLOAD_FOLDER"] = os.path.join("static", "icons")

Compress(app)

lock = Lock()
task_queue = queue.Queue()
task_status = {"specs": []}

from app.task_worker import task_worker

worker_thread = Thread(
  target=task_worker,
  daemon=True)
worker_thread.start()

from app.routes import (
  assoc,
  error,
  index,
  am_res,
  am_opts,
  assoc_qc,
  assoc_res,
  assoc_filt,
  assoc_opts,
  assoc_locus)

if __name__ == "__main__":
    app.run(
      debug=True,
      port=8888)#5101)

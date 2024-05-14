#!/usr/bin/env python

from flask import Flask
from flask_compress import Compress
from flask_executor import Executor

from app.dirs import (
  root_dir,
  data_dir,
  tmp_dir)

app = Flask(__name__)
app.secret_key = "275b250b2816223102a047d297676ce0"
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024

Compress(app)
executor = Executor(app)

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
      port=8888)

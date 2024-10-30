#!/usr/bin/env python

import os
from time import time

from app.dirs import tmp_dir


def clean():
    if os.path.exists(tmp_dir):
        for fn in os.listdir(tmp_dir):
            # 1 day = 86,400 seconds
            fn_path = os.path.join(
              tmp_dir,
              fn)
            if time() - os.path.getmtime(fn_path) > 86400:
                os.remove(fn_path)
    else:
        os.mkdir(tmp_dir)
    if os.path.exists("output"):
        for fn in os.listdir("output"):
            # 1 day = 86,400 seconds
            fn_path = os.path.join(
              "output",
              fn)
            if time() - os.path.getmtime(fn_path) > 86400:
                os.remove(fn_path)

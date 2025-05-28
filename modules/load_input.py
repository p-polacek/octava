#!/usr/bin/env python

import os
import json
import portalocker
import configparser

from app.dirs import data_dir

def load_gdata(genome):
    gdata_path = os.path.join(
      data_dir,
      genome,
      ".gdata")
    with open(gdata_path) as locked_f:
        portalocker.lock(locked_f, portalocker.LOCK_SH)
        contents = locked_f.read()
        portalocker.unlock(locked_f)
    gdata = json.loads(contents)
    return gdata

def load_config(genome):
    config = configparser.ConfigParser()
    conf_path = os.path.join(
      data_dir,
      genome,
      "config.ini")
    with open(conf_path) as locked_f:
        portalocker.lock(locked_f, portalocker.LOCK_SH)
        _ = config.read_file(locked_f)
        portalocker.unlock(locked_f)
    return config["setup"]

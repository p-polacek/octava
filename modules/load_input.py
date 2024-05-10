#!/usr/bin/env python

import os
import json
import configparser

from app.dirs import data_dir

def load_gdata(genome):
    gdata_path = os.path.join(
      data_dir,
      genome,
      ".gdata")
    with open(gdata_path, "r") as gdata_in:
        gdata = json.load(gdata_in)
    return gdata

def load_config(genome):
    config = configparser.ConfigParser()
    conf_path = os.path.join(
      data_dir,
      genome,
      "config.ini")
    _ = config.read(conf_path)
    return config["Configuration"]

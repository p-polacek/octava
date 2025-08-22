#!/usr/bin/env python

import os
import json
import portalocker
import configparser

from app.dirs import data_dir, tmp_dir


def load_gdata(genome, load_ext=False):
    if load_ext:
        fn_data_dir = tmp_dir
    else:
        fn_data_dir = data_dir

    gdata_path = os.path.join(fn_data_dir, genome, "stats", ".gdata")
    with portalocker.Lock(gdata_path, "r", flags=portalocker.LOCK_SH) as locked_f:
        contents = locked_f.read()
    gdata = json.loads(contents)
    return gdata


def load_config(genome, load_ext=False):
    if load_ext:
        fn_data_dir = tmp_dir
    else:
        fn_data_dir = data_dir

    config = configparser.ConfigParser()
    conf_path = os.path.join(fn_data_dir, genome, "stats", "config.ini")
    with portalocker.Lock(conf_path, "r", flags=portalocker.LOCK_SH) as locked_f:
        _ = config.read_file(locked_f)
    return config["setup"]

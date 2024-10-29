#!/usr/bin/env python

import os
import sys

curr_dir = os.path.dirname(
  os.path.realpath(
    __file__))
root_dir = os.path.abspath(
  os.path.join(
    curr_dir,
    ".."))
sys.path.insert(
  0,
  root_dir)

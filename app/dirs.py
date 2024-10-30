#!/usr/bin/env python
import os

root_dir = ""
data_dir = "data"
tmp_dir = "tmp"

tool_dir = os.path.join(
  root_dir,
  "modules",
  "tools")


plink1_path = os.path.join(
  tool_dir,
  "plink")
plink2_path = os.path.join(
  tool_dir,
  "plink2")
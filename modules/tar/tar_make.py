#!/usr/bin/env python

import os
import sys
import json
import tarfile

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from app.dirs import tmp_dir, data_dir

def create_tar(files, task_id):
    print("<b>Preparing input files</b>")
    tar_path = os.path.join(data_dir, f"{task_id}.in.tar.gz")
    with tarfile.open(tar_path, "w:gz") as tar:
        for rel_path in files:
            abs_path = os.path.join(data_dir, rel_path)
            if os.path.exists(abs_path):
                tar.add(abs_path, arcname=rel_path)
            else:
                raise OSError(f"Missing file: {rel_path}")
        in_params = f"{task_id}.cmd"
        tar.add(os.path.join(tmp_dir, in_params), arcname=in_params)
    print(f"<b>Your files are ready at {tar_path}</b>")


if __name__ == "__main__":
    files, task_id = sys.argv[1:3]
    files = json.loads(files)
    create_tar(files, task_id)

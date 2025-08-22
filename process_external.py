#!/usr/bin/env python

import os
import sys
import json
import shutil
import signal
import tarfile
import subprocess

from app.dirs import tmp_dir, data_dir

def process_ext(mem):
    ext_files = [
        f for f in os.listdir(data_dir)
        if f.endswith(".in.tar.gz")]

    for tar_file in ext_files:
        print(f"Extracting {tar_file}")
        cmd_file = tar_file.replace(".in.tar.gz", ".cmd")
        tar_path = os.path.join(data_dir, tar_file)
        with tarfile.open(tar_path, "r:gz") as tar:
            names = [m.name.split("/")[0] for m in tar.getmembers() if m.name]
            top_name = list(dict.fromkeys(names))[0]

            tar.extractall(path=tmp_dir)

        with open(os.path.join(tmp_dir, cmd_file)) as cmd_in:
            cmd = json.load(cmd_in)
            cmd.append(mem)
        pr = subprocess.Popen(
            cmd,
            preexec_fn=os.setsid,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        for line in pr.stdout:
            print(line, end="")
        try:
            pr.wait()
        finally:
            try:
                if pr.poll() is None:
                    os.killpg(os.getpgid(pr.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass

    shutil.rmtree(os.path.join(tmp_dir, top_name))

if __name__ == "__main__":
    mem = sys.argv[1] if len(sys.argv) == 2 else 0
    process_ext(mem)
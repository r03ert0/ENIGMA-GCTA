#!/usr/bin/env python3

import subprocess
import os

def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()


config_file = os.path.join(os.path.dirname(__file__), "config_test.yml")
src_dir = os.path.join(getannexdir(), "src")

def test_pipeline():
    subprocess.run([os.path.join(src_dir, "00.run_all.py"), config_file], check=True)
     

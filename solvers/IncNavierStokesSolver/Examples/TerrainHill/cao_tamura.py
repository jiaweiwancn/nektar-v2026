#!/usr/bin/env python3
"""Public entry point for the versioned Cao--Tamura experimental workflow."""
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'Tests/RANS'))
from cao_tamura_experiment import main

if __name__ == '__main__':
    main()

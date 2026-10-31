#!/usr/bin/env python3
"""Prepare a terrain-wall session using the versioned verification generator."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'Tests/RANS'))
from generate_terrain import main
if __name__ == '__main__':
    main()

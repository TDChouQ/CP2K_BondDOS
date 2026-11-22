"""Entry point: no arguments -> GUI, otherwise command line (see `python main.py -h`)."""
import sys
import multiprocessing

from bonddos.cli import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main())

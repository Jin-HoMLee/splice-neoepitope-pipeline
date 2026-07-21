import os
import sys

# Make the harness package importable when pytest is pointed at this directory.
sys.path.insert(0, os.path.dirname(__file__))

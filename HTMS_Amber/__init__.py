import os
import sys

current_package_dir = os.path.dirname(os.path.abspath(__file__))

if current_package_dir not in sys.path:
    sys.path.insert(0, current_package_dir)

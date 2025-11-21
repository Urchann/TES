#!/usr/bin/env python

r'''
 _    _      _        _                                                                                                           _                    
| |_ | |__  (_)  ___ | | __ _ __    ___  ___  ___         _ __ ___    ___   __ _  ___  _   _  _ __   ___  _ __ ___    ___  _ __  | |_     _ __   _   _ 
| __|| '_ \ | | / __|| |/ /| '_ \  / _ \/ __|/ __|       | '_ ` _ \  / _ \ / _` |/ __|| | | || '__| / _ \| '_ ` _ \  / _ \| '_ \ | __|   | '_ \ | | | |
| |_ | | | || || (__ |   < | | | ||  __/\__ \\__ \       | | | | | ||  __/| (_| |\__ \| |_| || |   |  __/| | | | | ||  __/| | | || |_  _ | |_) || |_| |
 \__||_| |_||_| \___||_|\_\|_| |_| \___||___/|___/ _____ |_| |_| |_| \___| \__,_||___/ \__,_||_|    \___||_| |_| |_| \___||_| |_| \__|(_)| .__/  \__, |
                                                  |_____|                                                                                |_|     |___/ 
'''

__version__ = '1.1' # 2025.12.31
__author__  = 'Ryota Fukuda' # mailto:25la018c@rikkyo.ac.jp

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
import pandas as pd
import argparse
import os
import subprocess
import shutil
from scipy.ndimage import uniform_filter1d
import re

class Thickness:
    def __init__(self, dattxt):
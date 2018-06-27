#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:25:36 2018

@author: brunella
"""

from __future__ import absolute_import

import sys
sys.path.append('/Users/brunella/Desktop/Simulazioni_PythonProper/MetisCoronagraphSimulator/CATS/cats/cats_simus')
from multi_prop import multi_prop

pixelsize = 5.
beam_ratio = pixelsize*4.85e-9/(3.8*1e-6/diam)

multi_prop(1024,3.8, 'prova1', RAVC=True, LS=True)



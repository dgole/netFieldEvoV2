#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
from matplotlib.backends.backend_pdf import PdfPages
import time
################################################################################
fig = plt.figure(figsize=(6,4), dpi=100)
ax = []
ax.append(plt.subplot2grid((1, 1), (0, 0), rowspan=1))
ax[0].plot(x, y, color=(0,0,1,1), linewidth=1)
################################################################################








################################################################################
fig.tight_layout()
fig.savefig('./diagram_horizontal.png'); plt.clf()
################################################################################

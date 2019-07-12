#!/usr/bin/python
import numpy as np
import os
import sys
import resource
import time


a = np.asarray([0,1,2,3,4,5,6,7,8,9])
print(a)
print(a<8)
print(a>1)
print(a==5)
testArray = np.logical_and(a<8, a==5)
print(testArray*3.0)

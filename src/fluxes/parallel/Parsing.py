import numpy as np
import math
import pylab as py
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma
from decimal import Decimal
import argparse

parser = argparse.ArgumentParser(description='Choose between two functions.')

parser.add_argument('function', help='function name. Options to use: first(a,b) , second(a,b), third, fourth, fifth, sixth, seventh, eighth, ninth, tenth, eleventh, twelveth...')
parser.add_argument('x', type=float,help='first argument passed')
parser.add_argument('y', type=float,help='second argument passed')
args = parser.parse_args()

def first(a,b):
    answer = a+b
    print("first function, sum = ", answer)

def second(a,b):
    answer = a*b
    print("second function, product = ",answer)

locals()[args.function](args.x, args.y)
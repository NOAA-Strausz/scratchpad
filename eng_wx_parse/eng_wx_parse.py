#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 14:21:14 2019

@author: dave
"""

filename = 'CM002_02_04_2019'

file = open(filename).read()

test = file.split('%%')
wind=[]
atrh=[]
baro=[]
for x in test:
    lines = x.split('\n')
    if lines[0] == 'WIND':
        if wind == []:
            wind=wind+lines[1:-1]
        else:
            wind=wind+lines[2:-1]
    if lines[0] == 'ATRH':
        if atrh == []:
            atrh = atrh + lines[1:-1]
        else:
            atrh = atrh + lines[2:-1]
    if lines[0] == 'BARO':
        if baro == []:
            baro = baro + lines[1:-1]
        else:
            baro = baro + lines[2:-1]

for x in wind:
    print(x)

for x in baro:
    print(x)

for x in atrh:
    print(x)
            
    

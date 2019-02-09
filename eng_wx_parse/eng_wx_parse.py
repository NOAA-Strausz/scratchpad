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
for x in test:
    lines = x.split('\n')
    if lines[0] == 'WIND':
        if wind == []:

            wind=wind+lines[1:-1]
        else:
            wind=wind+lines[2:-1]


for x in wind:
    print(x)


            
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:25:43 2019

@author: strausz
"""

import math
from haversine import haversine



def hav(point1, point2):
    
    lat1, lon1 = point1
    lat2, lon2 = point2
    
    r = 6371 #radius of earth in km
    
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    d = 2*r*math.asin(math.sqrt((math.sin(dlat/2))**2 + math.cos(lat1)*math.cos(lat2)*(math.sin(dlon/2))**2))
    return d

def destination(lat1, lon1, brg, dist, nm):
    
    if nm:
        r=3440
    else:
        r=6371
    
    lat2 = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                     math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                     math.cos(math.radians(brg)))
    
    lon2 = math.radians(lon1) + math.atan2(math.sin(math.radians(brg)) *
                        math.sin(dist/r) * math.cos(math.radians(lat1)),
                        math.cos(dist/r) - math.sin(math.radians(lat1)) *
                        math.sin(math.radians(lat1)))
    
    return(math.degrees(lat2), math.degrees(lon2))
    

def find_box(lat1, lon1, dist, nm):
    
    if nm:
        r=3440
    else:
        r=6371
    
    dist = dist/2
    
    wlon = math.radians(lon1) + math.atan2(math.sin(math.radians(270)) *
                        math.sin(dist/r) * math.cos(math.radians(lat1)),
                        math.cos(dist/r) - math.sin(math.radians(lat1)) *
                        math.sin(math.radians(lat1)))
    elon = math.radians(lon1) + math.atan2(math.sin(math.radians(90)) *
                        math.sin(dist/r) * math.cos(math.radians(lat1)),
                        math.cos(dist/r) - math.sin(math.radians(lat1)) *
                        math.sin(math.radians(lat1)))
    nlat = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                     math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                     math.cos(math.radians(0)))
    slat = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                     math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                     math.cos(math.radians(180)))
    
    wlon = round(math.degrees(wlon), 3)
    elon = round(math.degrees(elon), 3)
    nlat = round(math.degrees(nlat), 3)
    slat = round(math.degrees(slat), 3)
    
    return([nlat, wlon], [nlat, elon], [slat, wlon], [slat, elon])
    
    
    
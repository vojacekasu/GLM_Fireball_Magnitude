# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 10:09:48 2022

@author: Vlastimil Vojacek
"""

import matplotlib.pyplot as plt
import pandas as pd 
import math
#from math import sin, cos, sqrt, atan2, radians, pi
import pyproj
import re
#import numpy as np

###############    provide these information before run:
csvfile = "GLM_csv_file.csv"  #insert the CSV file containig GLM data
height = 16000  #insert height of fireball in atmosphere in meters, if provided
v = 20      #insert velocity of the fireball in km/s
plt.ylim([-14, -26])  #set the magnitude range of the plot if needed
############################################################################

head = []
#reads header of the csv file
with open(csvfile, 'r') as f:
    for line in f:
        if line.startswith('#'):
            head.append(line[1:].strip().split(','))
        else:
            break #stop when there are no more #

#takes latitude, longitude, distance to Earth of the GLM satellite
SatInfo = head[4][1].split('/')
SatLatitude = float(SatInfo[0])
SatLongitude = float(SatInfo[1] )
SatDistance = SatInfo[2]
SatDistance = re.findall('\d*\.?\d+',SatDistance)
SatDistance = float(SatDistance[0])*1000  #converts distance from km to m



#reads rest of the csv file as dataframe, first row are names of columns
df = pd.read_csv(csvfile, delimiter=",", comment='#')

#to start the time at zero
start = df['time (ms)'][0]
df['time (s)'] = (df['time (ms)'])/1000

#converts time from GLM in miliseconds into date time of the event
#date0 = datetime.datetime(1970, 1, 1, 0,0,0)
#delta = timedelta(seconds = start/1000)
#eventTime =  date0 + delta



#time = eventTime.time()
#seconds = time.second #takes only seconds from date and time of the event

df['datetime'] = pd.to_datetime(df["time (ms)"], unit='ms') #converts time in GLM in sec to date time of event (miliseconds from 1970 01 01 00h 00m 00s)
df['seconds'] = df['datetime'].dt.second +  df['datetime'].dt.microsecond/1000000  #uses only seconds from the date time of the event

############################################################################


#############from longitude, latitude computes distance from fireball to satellite
def distanceToGLM (lon,lat,hae):
#####compute  GLM position from lon, lat, high above sea level (hae) to x,y,z coordinates  
    latGLM = SatLatitude
    lonGLM = SatLongitude
    haeGLM = SatDistance
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    xGLM ,yGLM, zGLM = transformer.transform(lonGLM,latGLM,haeGLM,radians = False)

#####computes  meteor position from lon, lat, high above sea level (hae) to x,y,z coordinates  

    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    x ,y, z = transformer.transform(lon,lat,hae,radians = False)

    
 ###computes distance in km between meteor and GLM   
    dist = (math.sqrt( (x - xGLM)**2 + (y - yGLM)**2  + (z - zGLM)**2))/1000
    return dist


df['meteor_to_GLM_dist_(km)'] = df.apply(
    lambda row: distanceToGLM (
        lon = row['longitude'],
        lat = row['latitude'],
        hae = height
    ),
    axis=1
)
############################################################################


############computes magnitude from the GLM reported energy


dt = 0.002               #frequency of the detector seconds
A = math.pi*0.0000558**2 #effective lens aperture

mag = []

#parameters from the dependence of the signal at 777nm on the meteor velocity
a = 0.0948
c = -3.5

def magnitude (velocity, energ, dist, dt):     #computes absolute magnitude
    m = a * velocity - 2.5*math.log10(energ*dist**2) + 2.5*math.log10(dt*A) + c 
    return m



df['mag'] = df.apply(
    lambda row: magnitude (
        velocity = v,
        energ = row['energy (joules)'],
        dist = row['meteor_to_GLM_dist_(km)'],
        dt = dt
    ),
    axis=1
)
    
    

############################################################################

##### plots the lightcurve

x = df['seconds'] 
y = df['mag']

LC = plt.plot(x,y, label='GLM lightcurve in mag', color = 'red', linestyle='-', marker='o')


######correction for the brightness of the meteor

m0 = -11          #'zero' brightness in mag from which the correction is computed

plt.tick_params(labelsize=30)
plt.xlabel("time [s]", size = 30)
plt.ylabel("mag", size = 30)
plt.legend(loc='upper left',prop={"size":30})

############################################################################




df.to_csv('output.csv')  ##writes a file with computed values
plt.show()
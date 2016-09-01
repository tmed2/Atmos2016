# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:11:35 2016

@author: tmed2

Simple script to take Nasa satellite data and upack it into a binary file
in a more regular format. The data in question come from the joined data (in
the .csv files) from:
http://acdb-ext.gsfc.nasa.gov/Data_services/cloud_slice/new_data.html
"""


import csv
import numpy as np


gas_conc = np.array([])
longitude = np.array([])
latitude = np.array([])

with open('NASA/NasaSatellite/tco_oct04_to_dec13.csv', 'r') as file:
    csvfile = csv.reader(file)
    for index, row in enumerate(csvfile):
        #for reasons, the values are not separated by a fixed char or string
        row = row[0].split()
        conc = float(row[0])
        long = float(row[1])
        lat = float(row[2])
        gas_conc = np.append(gas_conc, conc)
        longitude = np.append(longitude, long)
        latitude = np.append(latitude, lat)


dates = []
start_month = 10
for month_since in range(121):
    #probably could find a more elegant way
    d = 15
    m =  (start_month + month_since) % 12
    y = 2004 + ((start_month + month_since) // 12)
    if m == 0:
        m = 12
        y -= 1
    dates.append(str(d) + '/' + str(m) + '/' + str(y))
    
    
#Assuming the .csv data are ordered...
with open('NASA/NasaSatellite/NasaDataDU.dat', 'ab+') as file:
    #lets hope there are no missing values    
    repeat = 2592
    for index, conc in enumerate(gas_conc):
        d_index = index // repeat
        row = dates[d_index]  + ' ' + str(conc) + ' ' + str(latitude[index])
        row = row + ' ' + str(longitude[index]) + '\n'
        file.write(row.encode())
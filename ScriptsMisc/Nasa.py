# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:39:35 2016

@author: tmed2

NB this is deprecated as of 29th July 2016. Use the AtmosGasTools package.

This script contains functions which search the Nasa data (which have been 
pre-processed with 'NasaTropospheric.py') for a given position. N.B. Double
check the file paths before use!
"""


import numpy as np
import datetime
import AtmosGasTSclass as ag


def nasa_reader(latitude, longitude):
    """
    Searches the NasaData file (assumed to be in "Data/Nasa/NasaData.dat") and
    returns a time series of data for the location at 'latitude' and 
    'longitude'. Will not return negative concentration values.
    
    latitude -- float; latitude of the site of interest, North is positive
                and South negative.
                
    longitude -- float; longitude of the site of interest, East is positive,
                 and West is negative.
    """
    
    conc_arr = np.array([])
    date_arr = np.array([])
    with open('NASA/NasaSatellite/NasaDataDU.dat', 'rb') as file:
        for row in file:
            str_row = row.decode()
            list_row = str_row.split()
                
            date = list_row[0]
            conc = float(list_row[1])
            lat = float(list_row[2])
            long = float(list_row[3])
            
            if (lat == latitude) and (long == longitude):
                date = date.split('/')
                d = int(date[0])
                m = int(date[1])
                y = int(date[2])
                
                datetime_obj = datetime.datetime(y, m, d)
                
                if conc >= 0:
                    conc_arr = np.append(conc_arr, conc)
                    date_arr = np.append(date_arr, datetime_obj)
                else:
                    pass
    
    if (conc_arr.size == 0) and (date_arr.size == 0):
        print('There is no data, check latitude and longitude')
        return None
    
    return [date_arr, conc_arr]

def near_loc(i_latitude, i_longitude):
    """
    returns the values of the lattiude and longitude, with extant data, closest
    to i_lattitude and i_longitude. If the chosen point is equidistant to more 
    than one point in the data then the first latitude and longitude with this 
    property will be returned. NB, no periodic correction is applied!
    
    i_lattitude -- float; latitude of the site of interest, North is positive
                   and South negative.
                   
    i_longitude -- float; longitude of the site of interest, East is positive,
                   and West is negative.
    """
    lat_arr = np.array([])
    long_arr = np.array([])
    
    date_start = '15/10/2004'
    date = date_start #?
    
    #this depends on the data being organised into blocks coressponding to
    #a single date, or at least at for the first date (date_start)
    with open('NASA/NasaSatellite/NasaDataDU.dat', 'rb') as file:    
        for row in file:
            str_row = row.decode()
            list_row = str_row.split()
                    
            date = list_row[0]
            #this way each position is only counted once
            if date != date_start:
                break
            else:
                pass
                
            lat = float(list_row[2])
            long = float(list_row[3])
                
            lat_arr = np.append(lat_arr, lat)
            long_arr = np.append(long_arr, long)
    
    dist_lat = lat_arr - i_latitude
    dist_long = long_arr - i_longitude
    dist_lat = np.abs(dist_lat)
    dist_long = np.abs(dist_long)
    
    latitude = lat_arr[dist_lat.argmin()]
    longitude = long_arr[dist_long.argmin()]
    return latitude, longitude

"""
#Here is an example of use:
These are the coordinates of Heimaey, Iceland
lat, lng = near_loc(63.4277, -20.2674)
data = nasa_reader(lat, lng)

h = ag.AtmosGasTS('O3', 'ppb', 'Heimaey Nasa Ozone', 
                  'Heimaey Nasa Ozone Spectrum', data[0], data[1], 'Nasa')
h.fast_find_spectrum(25000, -8.5, -6.6847)

h.ts_plot(join_data = True, save = True)
h.ps_plot(join_data = True, days = True, save = True, rel = False)
h.ps_plot(join_data = True, days = True, save = True, rel = True)
"""
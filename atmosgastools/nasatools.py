# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:39:35 2016

@author: tmed2

This module contains functions designed to search and process tropospheric
column ozone data from Nasa satellites:
http://acdb-ext.gsfc.nasa.gov/Data_services/cloud_slice/new_data.html

NB:
#File paths must be relative to the calling script (not this, the defining
 script)


Available functions:
-unpack_nasa_csv: decodes .csv data into an expicit time series in a .dat file
-nasa_reader: reads in the time series of a particular location
-near_loc: given a location, finds the nearest location with a time series
"""


import csv
import datetime

import numpy as np




#this is a very shoddy wrapper of a quick script, generalise later...
def unpack_nasa_csv():
    """
    take Nasa satellite data and upack it into a binary file in a more regular
    format. The data in question come from the joined data (in the .csv files)
    from:
    http://acdb-ext.gsfc.nasa.gov/Data_services/cloud_slice/new_data.html
    """
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
        m = (start_month + month_since) % 12
        y = 2004 + ((start_month + month_since) // 12)
        if m == 0:
            m = 12
            y -= 1
        dates.append(str(d) + '/' + str(m) + '/' + str(y))


    #Assuming the .csv data are ordered...
    with open('NASA/NasaSatellite/NasaDataDU.dat', 'ab') as file:
        #lets hope there are no missing values
        repeat = 2592
        for index, conc in enumerate(gas_conc):
            d_index = index // repeat
            row = dates[d_index]  + ' ' + str(conc) + ' ' + str(latitude[index])
            row = row + ' ' + str(longitude[index]) + '\n'
            file.write(row.encode())
    return None


def nasa_reader(latitude, longitude):
    """
    Searches the NasaData file (assumed to be in "Data/Nasa/NasaData.dat") and
    returns a time series of data for the location at 'latitude' and
    'longitude'. Will not return negative concentration values. NB this given
    position must be in the data set, if it is not find the nearest position
    with 'near_loc()'

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
    date = date_start  # ?

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
import AtmosGasTools.AtmosGasTS as ag
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
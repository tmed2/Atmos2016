# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 11:54:48 2016

@author: tmed2

NB this may throw errors as it used an older version of AtmosGasTSclass;
minor changes may be needed to make it work
"""

import csv
import datetime
import numpy as np
#import matplotlib.pyplot as plt
#from scipy import signal
import AtmosGasTSclass as ag


def lhr_reader(filename, species_col, baseline = False):
    """
    Reads in the .csv file. The files must be in the path: Data/HeathrowAirport
    in the script directory, which must be formated like the Heathrow airport 
    data. Returns a time series of the concentrations.
    
    Arguments:
    filename -- string; name of the .csv file (without the file ending) in the
                "Data" folder.
    
    species_cols -- array, of two ints; the column index of the .csv file which 
                    contains the relevant concentration data, and its baseline
                    respectively. If Baseline is false, a 1 element array
                    should work also.
    
    Optional Arguments:
    baseline -- Boolean, default is False; if true, the 'baseline' of the gas 
                is subtracted from the usual data. Otherwise raw concentration 
                data are used.
    """
    
    time_array = np.array([])
    conc_array = np.array([])
    with open('Data/HeathrowAirport/' + filename + '.csv', 'r') as file:
        data = csv.reader(file, delimiter=',')
        for index, row in enumerate(data):
            if row[0] == 'date':
                print('No data; line -', index + 1)
            else:
                #checks the baseline subtraction, may throw errors if the data
                #contains strings below the first line
                if baseline:
                    conc = float(row[species_col[0]]) - float(row[species_col[1]])
                else:
                    conc = float(row[species_col[0]])
                
                #ensures only, positive definite data
                if (conc >= 0):
                    
                    date_and_time = row[0].split(' ')
                    
                    date = date_and_time[0].split('/')
                    time = date_and_time[1].split(':')
                    
                    year = int(date[2])
                    month = int(date[1])
                    day = int(date[0])
                    
                    if int(time[0]) == 24:
                        hour = 0
                    else:
                        hour = int(time[0])
                    minute = int(time[1])
                    second = int(time[2])
                    
                    datetime_obj = datetime.datetime(year, month, day, hour, 
                                                     minute, second)
                    time_array = np.append(time_array, datetime_obj)
                    conc_array = np.append(conc_array, conc)
                else:
                    print('Negative concentration, not included; line -', 
                          index + 1)
                    pass
    return [time_array, conc_array]

start_comp = datetime.datetime.now()
co_data= lhr_reader('SITE32_bl_local_raw_041012_121112', [3], 
                    baseline = False)

co = ag.AtmosGasTS('CO', 'ppb', 'Heathrow32 CO Time Series 2 1hr Averaging', 
                   'Heathrow32 CO Frequency Spectrum 2 1hr Averaging', 
                   co_data[0], co_data[1], subpath = 'HeathrowAirport')

#in this data set, each step is 20s
co.moving_mean(180)          
co.ts_plot(join_data = False, save = True)
co.fast_find_spectrum(50000, -6.414, -3.8)
co.ps_plot(join_data = False, save = True, days = True)
end_comp = datetime.datetime.now()
print('Time elapsed:', end_comp - start_comp)
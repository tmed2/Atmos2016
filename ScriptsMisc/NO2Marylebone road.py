# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:34:45 2016

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

def no2_reader(filename):
    """
    Reads in the .csv file. The files must be in a folder called "Data" in
    the script directory.
    """
    
    time_array = np.array([])
    conc_array = np.array([])
    with open('Data/' + filename + '.csv', 'r') as file:
        data = csv.reader(file, delimiter=',')
        for index, row in enumerate(data):
            
            if row[2] == 'No data':
                print('No data; line -', index)
            #ensures only verified, positive definite data
            elif (float(row[2]) >= 0) and (row[3][0] == 'V'):
                conc_array = np.append(conc_array, row[2])
            
                date = row[0].split('-')
                time = row[1].split(':')
                
                year = int(date[0])
                month = int(date[1])
                day = int(date[2])
                
                if int(time[0]) == 24:
                    hour = 0
                else:
                    hour = int(time[0])
                minute = int(time[1])
                second = int(time[2])
                
                datetime_obj = datetime.datetime(year, month, day, hour, 
                                                 minute, second)
                time_array = np.append(time_array, datetime_obj)
            else:
                pass
    return [time_array, conc_array]


no2data = no2_reader('NO2trimmed')
no2 = ag.AtmosGasTS('NO2', 'ugm^-3', 'NO Time Series',
                    'NO Frequency Spectrum', no2data[0], no2data[1])
no2.find_spectrum(20000, -7.5, -4.3)
no2.spec_plot(join_data = True, days = True)
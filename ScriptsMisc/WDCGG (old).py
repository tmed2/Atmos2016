# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:25:14 2016

@author: tmed2

The purpose of this script is to take .dat files downloaded from the WDCGG
(http://ds.data.jma.go.jp/gmd/wdcgg/wdcgg.html) and parse them for relevant 
data. If time permits, URL requests may be implemented in order to automate
the data gathering process.

A range of plotting tools, and some statistical functions (e.g. a power 
spectrum estimator) have also been implemented.
"""


import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def dat_reader(filename):
    """
    Reads in the .dat file. The .dat files must be in a folder called "Data" in
    the script directory.
    """
    
    header = []
    data = []
    with open('Data/' + filename + '.dat', 'rb') as file:
        for row in file:
            string_row = row.decode()
            if string_row[0] == 'C':
                header.append(string_row)
            else:
                data.append(string_row.split())
            
    return [header, data]
    

def data_joiner(data0, data1):
    """
    Takes two data matricies, like those produced by dat_reader, and
    returns a new matrix consisting of the rows od 'data1' appended to
    those of .data0'.
    """
    joined_data = data0
    for row in data1:
        joined_data.append(row)
    return joined_data


def data_md(data):
    """
    Takes in a 'data' matrix produced by dat_reader, and returns a more
    plotable format. This function is desgined for data recorded on a monthly
    or daily basis.
    """
    
    gas_conc = np.array([])
    date_array = np.array([])
    
    for index, row in enumerate(data):
        #care, there are two date columns that may be swtiched
        date = row[0].split('-')
        
        year = int(date[0])
        month = int(date[1])
        day = int(date[2])
        
        date_obj = datetime.date(year, month, day)
        gas = float(row[4])
        
        if gas >= 0:
            gas_conc = np.append(gas_conc, gas)
            date_array = np.append(date_array, date_obj)
        elif gas < 0:
            #line starts at the date in the .dat file, it does not include the
            #header
            print('Gas concentration is negative; not included. line -', index + 1)
        else:
            #line starts at the date in the .dat file, it does not include the
            #header
            print('There is some problem; not included. line -', index + 1)
        
        
    time_series = np.array([date_array, gas_conc])
    return time_series


#this should also work with monthly and daily data
def data_h(data):
    """
    Does the same thing as data_md, but uses datetime objects (instead of
    date objects) in the time series
    """
    
    gas_conc = np.array([])
    datetime_array = np.array([])
    
    for index, row in enumerate(data):
        #care, there are two date and time columns that may be swtiched
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
        
        datetime_obj = datetime.datetime(year, month, day, hour, minute)
        gas = float(row[4])
        
        if gas >= 0:
            gas_conc = np.append(gas_conc, gas)
            datetime_array = np.append(datetime_array, datetime_obj)
        elif gas < 0:
            #line starts at the date in the .dat file, it does not include the
            #header
            print('Gas concentration is negative; not included. line -', index + 1)
        else:
            #line starts at the date in the .dat file, it does not include the
            #header
            print('There is some problem; not included. line -', index + 1)
        
        
    time_series = np.array([datetime_array, gas_conc])
    return time_series


def deviation_check(time_series, sigma_cut):
    """
    Returns a time series with values that are 'sigma_cut' standard 
    deviations from the mean removed. 'time_series' must be like those produced
    by the 'data_monthly' function, ie rows of [time, gas_conc].
    """
    
    average_gas_conc = np.mean(time_series[1])
    sigma_gas_conc = np.std(time_series[1])
    upper_cut = average_gas_conc + sigma_cut*sigma_gas_conc
    lower_cut = average_gas_conc - sigma_cut*sigma_gas_conc
    
    new_gas_conc = np.array([])
    new_time = np.array([])
    
    for index, conc in enumerate(time_series[1]):
        if (conc < upper_cut) and (conc > lower_cut):
            new_gas_conc = np.append(new_gas_conc, conc)
            new_time = np.append(new_time, time_series[0][index])
        elif (conc > upper_cut) or (conc < lower_cut):
            #this concerns the input time series, index starts at 0
            print('Index ', index, ' removed; concentration out ofrange')
        else:
            #this concerns the input time series, index starts at 0
            print('Index ', index, ' removed; there is some problem')
    
    new_time_series = np.array([new_time, new_gas_conc])
    return new_time_series


def gas_data_plotter(times, gas_concs, conc_unit, gas_name, title, 
                     join_data = True, save = True):
    """ 
    'times' and 'gas_concs' must be mutually ordered arrays of numeric 
    types. The remaining arguments must be strings, except the last two
    which are Booleans.
    """
    
    plt.figure()
    if join_data:
        plt.plot(times, gas_concs)
    else:
        plt.plot(times, gas_concs, ls = '', marker = '+', ms = 1)
        
    plt.xlabel('Date', fontsize = 18)
    plt.ylabel(gas_name + ' Concentration / ' + conc_unit, fontsize = 18)
    plt.title(title, fontsize = 22)
    if save:
        plt.savefig('Figures/' + title + '.pdf', format = 'pdf')
        plt.close()
        print('Figure saved as: ', title)
        
    return None


#this is meant to be a quick fix, replace with better procedure later
def gas_multi_plot(multi_times, multi_gas_concs, name_array, conc_unit, 
                   gas_name, title):
                       
    plt.figure()
    for index, name in enumerate(name_array):
        plt.plot(multi_times[index], multi_gas_concs[index], label = name)
    plt.xlabel('Date', fontsize = 18)
    plt.ylabel(gas_name + ' Concentration / ' + conc_unit, fontsize = 18)
    plt.title(title, fontsize = 22)
    plt.legend(loc = 'best')
    plt.savefig(title + '.pdf', format = 'pdf')
    plt.close()
    return None


def periodogram(time_series, f_samps, mag_low, mag_high):
    """
    Finds the Lomb-Scargle periodogram of a 'time_series' like those produced
    by the data_x functions. 'f_samp' number of frequency samples are used
    in the range [10^mag_low, 10^mag_high]. Returns a normalised spectrum in
    the original units
    """
    
    assert mag_low < mag_high
    times = np.array([])
    angfreqs = 2*np.pi*np.logspace(mag_low, mag_high, f_samps)
    #the date/datetime objects must be in increasing order
    t0 = time_series[0][0]
    for t in time_series[0]:
        delta_t = t - t0
        times = np.append(times, delta_t.total_seconds())
    
    concs = time_series[1].astype(np.float64)
    concs = concs - np.mean(concs)
    pgram = signal.lombscargle(times, concs, angfreqs)
    N = times.shape[0]
    norm_pgram = np.sqrt(4*(pgram/N))
    
    norm_pow_spectrum = np.array([angfreqs/(2*np.pi), norm_pgram])
    return norm_pow_spectrum


def plot_freq_domain(pow_spectrum, conc_unit, freq_unit, gas_name, title,
                     join_data = True, save = True):
    """
    Plots a spectrum of those produced by 'periodogram'. 
    
    Positional arguments:
    pow_spectrum -- double array; frequencies and corresponding amplitudes
                    like those produced by 'periodogram'
    conc_unit -- string; unit of gas concentration
    freq_unit -- string; unit of frequency
    gas_name -- string; name of gas
    title -- string; will be the title on the graph
    
    keyword arguments:
    join_data -- Boolean; if true, the data are connected with lines. Otherwise
                 only points are plotted with + markers
    save -- Boolean; if true, the figure is saved as a .pdf with 'title' as the
            filename. Otherwise the figure will be kept open in the viewer for
            inspection and must be saved manually. There must be a folder
            called 'Figures' in the script directory.
    """
    
    plt.figure()
    if join_data:
        plt.plot(pow_spectrum[0], pow_spectrum[1])
    else:
        plt.plot(pow_spectrum[0], pow_spectrum[1], ls = '', marker = '+', ms = 1)
        
    plt.xlabel('frequency / ' + freq_unit, fontsize = 18)
    plt.ylabel(gas_name + ' Concentration / ' + conc_unit, fontsize = 18)
    plt.title(title, fontsize = 22)
    if save:
        plt.savefig('Figures/' + title + '.pdf', format = 'pdf')
        plt.close()
        print('Figure saved as: ', title)
    return None




"""
Izana_m = dat_reader('IzanaMonthly')
Izana_ts_m = data_md(Izana_m[1])
Izana_ts_m = deviation_check(Izana_ts_m, 3)

gas_data_plotter(Izana_ts_m[0], Izana_ts_m[1], 'ppb', 'Ozone',
                 'Izana Monthly Data', join_data = True, save = True)

Izana_spectrum = periodogram(Izana_ts_m, 40000, -10, -6.65)
plot_freq_domain(Izana_spectrum, 'ppb', 'hz', 'Ozone', 
                 'Frequency Spectrum Izana', join_data = False, save = True)
"""
heimaey93_h = dat_reader('Heimaey1993h')
heimaey94_h = dat_reader('Heimaey1994h')

heimaey9394_h = data_joiner(heimaey93_h[1], heimaey94_h[1])
heimaey9394_h = data_h(heimaey9394_h)
gas_data_plotter(heimaey9394_h[0], heimaey9394_h[1], 'ppb', 'Ozone',
                 'Heimaey 1993-1994', join_data = False, save = True)


heimaey9394_spectrum = periodogram(heimaey9394_h, 20000, -10, -4.9)
plot_freq_domain(heimaey9394_spectrum, 'ppb', 'hz', 'Ozone',
                 'Frequency Spectrum Heimaey 1993-1994', join_data = False,
                 save = True)

"""          
Tudor_m = dat_reader('TudorHillMonthly')
Tudor_ts_m  = data_md(Tudor_m[1])
Tudor_ts_m = deviation_check(Tudor_ts_m, 3)
gas_data_plotter(Tudor_ts_m[0], Tudor_ts_m[1], 'ppb', 'Ozone', 'Tudor Monthly')

Tudor_d = dat_reader('TudorHillDaily')
Tudor_ts_d  = data_md(Tudor_d[1])
Tudor_ts_d = deviation_check(Tudor_ts_d, 3)
gas_data_plotter(Tudor_ts_d[0], Tudor_ts_d[1], 'ppb', 'Ozone', 'Tudor Daily')

Eskdalemuir_m = dat_reader('EskdalemuirMonthly')
Eskdalemuir_ts_m = data_md(Eskdalemuir_m[1])
Eskdalemuir_ts_m = deviation_check(Eskdalemuir_ts_m, 3)

Izana_m = dat_reader('IzanaMonthly')
Izana_ts_m = data_md(Izana_m[1])
Izana_ts_m = deviation_check(Izana_ts_m, 3)

names = ['Tudor Hill', 'Eskdalemui', 'Izana']
multi_times_ex = np.array([Tudor_ts_m[0], Eskdalemuir_ts_m[0], Izana_ts_m[0]])
multi_gas_conc_ex = np.array([Tudor_ts_m[1], Eskdalemuir_ts_m[1], Izana_ts_m[1]])
gas_multi_plot(multi_times_ex, multi_gas_conc_ex, names, 'ppb', 'Ozone', 'Ozone Monthly')
"""
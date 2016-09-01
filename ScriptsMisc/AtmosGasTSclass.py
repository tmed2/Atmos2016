# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:12:56 2016

@author: tmed2

NB this is deprecated as of 29th July 2016. Use the AtmosGasTools package.

This script contains general functions and class defintions which are either
used by different scripts, or have definitions which are too long to place in
"ordinary" scripts.

NB some of the 'get...' methods were changed (in name only) to 'find...'
methods in the correlation class. Make sure that scripts which call these
methods call the 'find...' methods, if an AttributeError is thrown.
"""


import datetime
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import linregress
from scipy.interpolate import interp1d

from gatspy.periodic import LombScargleFast


#This is used to define the start of time, so that relative times will be
#consistent
EPOCH = datetime.datetime(1975, 1, 1)

###############################################################################

class AtmosGasTS(object):
    """
    Class containing attributes and methods relating to the time series
    of atmospheric gases.
    """
    def __init__(self, gas_name, gas_units, t_title, f_title, times, concs,
                 station, source):
        """
        Initialises the AtmosGasClass.

        gas_name -- string; name of the gas species, e.g. O3, NO2, CH4...

        gas_units -- string; units of concentration.

        t_title -- string; name to be used on time series.

        f_title -- string; name to be used on power spectra.

        times -- array; an array of date/datetime objects ordered with concs.

        concs -- array; an array of the gas concentrations corresponding to
                 the 'times' array.

        station -- string; file path of saved figures:
                   Stations/station/gas_name/title.pdf
        """


        self.g_name = gas_name
        self.g_unit = gas_units
        self.t_title = t_title
        self.f_title = f_title

        self.source = source
        self.station = station
        self.path = source + '/' + station + '/' + gas_name
        #these will be used to name the files, they are the same as x_title
        #but with newlines ('\n') removed.
        self.t_name = t_title.replace('\n', '')
        self.f_name = f_title.replace('\n', '')

        #containers for the raw data
        self.t_raw = times
        self.c_raw = concs.astype('float64')

        #these will become different to their raw counterparts if correction
        #methods are called; they are containers for processed data
        self.t = self.t_raw
        self.c = self.c_raw

        #will have to be calculated
        self.spectrum = None
        self.rel_spectrum = None

        #linear regression results from scipy.stats.linregress
        self.t_rel = None #since start, in seconds
        self.line_stats = None
        self.line = None

        #this will become a function of relarive times!
        self.interpol = None

        #These positions are optional and are set with 'set_global_position()'
        self.latitude = None
        self.longitude = None
        self.altitude = None
        return None

    @staticmethod
    def __t_rel(time_arr):
        """
        Returns an array of relative times (since EPOCH). This method was meant
        for use within other methods, hence privacy.

        time_arr -- array; an array of datetime objects whose relative times
                    will be found.
        """

        rel_time = np.array([])
        t0 = EPOCH
        for t in time_arr:
            delta_t = t - t0
            rel_time = np.append(rel_time, delta_t.total_seconds())
        return rel_time

    def find_t_rel(self, raw = False):
        """
        Finds the relative times (since the first datetime object).

        raw -- Boolean; if true the raw data are used. Otherwise the processed
               data are.
        """
        if raw:
            abs_time = self.t_raw
        else:
            abs_time = self.t
        self.t_rel = self.__t_rel(abs_time)
        return None

    def set_global_position(self, lat, long, alt = 0):
        """
        Sets the values of the latitude and longitude for this object
        """
        self.latitude = float(lat)
        self.longitude = float(long)
        self.altitude = float(alt)
        return None

    def cut_series(self, time_arr, raw = False):
        """
        Cuts out values from the time series that are not in the range (i.e.
        between the largest and smallest values) of 'time_arr'. The resulting
        arrays are placed in the 't' and 'c' processed data attributes.

        time_arr -- array; the largest and smallest values of this array will
                    be used as bounds on the data

        raw -- Boolean; if true the raw data are used. Otherwise the processed
               data are.
        """

        if raw:
            t = self.t_raw
            gas_conc = self.c_raw
        else:
            t = self.t
            gas_conc = self.c

        low = np.min(time_arr)
        high = np.max(time_arr)
        self.t = t[(t > low) & (t < high)]
        self.c = gas_conc[(t > low) & (t < high)]
        return None

    def find_line_stats(self, raw = False):
        """
        Finds the linear regression data using scipy.stats.linregress.

        raw -- Boolean; if true the raw data are used. Otherwise the processed
               data are. Make sure this is consistent with other method calls.
        """
        self.find_t_rel(raw = raw)
        if raw:
            gas_conc = self.c_raw
        else:
            gas_conc = self.c

        self.line_stats = linregress(self.t_rel, gas_conc)
        return None

    def find_line(self, raw = False):
        """
        Finds arrays of line values. These should be ordered with the 't' or
        't_raw' object arrays.

        raw -- Boolean; if true the raw data are used. Otherwise the processed
               data are. NB, this is only important if 'line_data' has not been
               calculated. It will have no effect if 'line_data' has already
               been calculated.
        """

        if self.line_stats is None:
            self.find_line_stats(raw)
        else:
            pass

        gradient = self.line_stats[0]
        intercept = self.line_stats[1]

        self.line = intercept + self.t_rel * gradient
        return None

    def write_reg_data(self, title = 'reg_stats'):
        """
        Writes the regression data into a text file. NB, if the filename and
        path are not changed then the data will be appended onto the existing
        file.

        title -- string; name of the file.
        """
        print('Regression data saved as:')
        print(self.path + '/' + title + '.txt')
        with open(self.path + '/' + title + '.txt',  'a') as txtfile:
            txtfile.write('These data were measured in ' + self.station + '\n')
            txtfile.write('Gradient is per second\n')
            txtfile.write('The gas is ' + self.g_name + '\n')
            txtfile.write('Regression statistics for the time series trend\n')
            txtfile.write('(gradient, intercept, r value, p value,')
            txtfile.write(' std error in gradient)' + '\n')
            txtfile.write(str(self.line_stats) + '\n\n\n')
        return None

    def find_interpol(self, raw = False, order = 'linear'):
        """
        Uses a scipy.interpolate class to find the interpolation function:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
        NB self.interpol is a function!

        raw -- Boolean; if true the raw data are used. Otherwise the processed
               data are used.

        order -- string (or int); value of the 'kind' parameter:
                 Specifies the kind of interpolation as a string (‘linear’,
                 ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’ where
                 ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline
                 interpolation of first, second or third order) or as an
                 integer specifying the order of the spline interpolator to
                 use. Default is ‘linear’.
        """
        if raw:
            times = self.t_rel
            concs = self.c_raw
        else:
            times = self.t_rel
            concs = self.c

        self.interpol = interp1d(times, concs, kind = order)
        return None

    def set_interpol_data(self, time_arr, raw = False, order = 'linear'):
        """
        Uses the 'interpol' method-attribute to find the value of the function
        at the points in time_arr. This method overwrites the the processed
        data attributes.

        time_arr -- array; array of times (including datetime objects) that is
                    used to calculate the interpolated data.

        raw -- Boolean; if True the raw times are used.

        order -- string (or int); value of the 'kind' parameter:
                 Specifies the kind of interpolation as a string (‘linear’,
                 ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’ where
                 ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline
                 interpolation of first, second or third order) or as an
                 integer specifying the order of the spline interpolator to
                 use. Default is ‘linear’.
        """
        if self.t_rel is None:
            self.find_t_rel(raw = raw)
        else:
            pass

        if self.interpol is None:
            self.find_interpol(raw = raw, order = order)
        else:
            pass

        self.t = time_arr
        time_arr = self.__t_rel(time_arr)
        self.c = self.interpol(time_arr)
        return None

    def regularise(self, start_ts, end_ts, interval, raw = False,
                   order = 'linear'):
        """
        Creates a series of gapless data

        start_ts -- datetime object; start of the regularised time series

        end_ts -- datetime object; nominal end of the regularised time series;
                  the number of intervals will be rounded down if 'interval' is
                  not an exact divisor of the total time difference.

        interval -- float; time, in hours, between consecutive points

        raw -- Boolean; if True the raw times are used.

        order -- string (or int); value of the 'kind' parameter:
                 Specifies the kind of interpolation as a string (‘linear’,
                 ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’ where
                 ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline
                 interpolation of first, second or third order) or as an
                 integer specifying the order of the spline interpolator to
                 use. Default is ‘linear’.
        """
        assert end_ts > start_ts
        #converts to seconds
        interval_s = 60*60*interval
        delta_tot = (end_ts - start_ts).total_seconds()
        interval_number = int(delta_tot//interval_s)
        new_t = np.array([])
        for i in range(interval_number):
            t_i = start_ts + datetime.timedelta(seconds = ((i+1)*interval_s))
            new_t = np.append(new_t, t_i)

        self.set_interpol_data(new_t, raw = raw, order = order)
        return None

    def reset_processed_series(self):
        """
        Resets the values of self.t and self.c to self.t_raw and self.c_raw.
        """
        self.t = self.t_raw
        self.c = self.c_raw
        return None

    def deviation_check(self, sigma_cut):
        """
        Removes values from self.c which are more than 'sigma_cut' standard
        deviations from the mean.

        sigma_cut -- float; the number of +/- standard deviations to include.
        """

        mean_gas_conc = np.mean(self.c)
        sigma_gas_conc = np.std(self.c)
        norm_gas_conc = (self.c - mean_gas_conc)/sigma_gas_conc
        norm_gas_conc = np.absolute(norm_gas_conc)

        self.t = self.t[norm_gas_conc < sigma_cut]
        self.c = self.c[norm_gas_conc < sigma_cut]
        return None

    def moving_mean(self, steps):
        """
        Averages groups of points, in the processed data, to form a new time
        series. If the data are missing values and/or have irregular sampling
        times, then the number of steps will not necessarily correspond to a
        fixed period of time!

        steps -- int; number of steps to be averaged at a time. The relevant
                 arrays will be trimmed of their ends to get an exact division
                 of steps.
        """

        concs = self.c
        times = self.t

        avrg_conc = concs[:(concs.size // steps) * steps].reshape((-1, steps))
        avrg_conc = np.mean(avrg_conc, axis = 1)

        times = times[:(times.size // steps) * steps].reshape((-1, steps))
        avrg_time = np.array([])
        #slow, think of a better way; dtatime objects are awkward
        for index, time_block in enumerate(times):

                #NB the times in between are not weighted into the mean/mid-
                #point time.
                start = time_block[0]
                end = time_block[-1]
                middle = start + (end - start)/2
                avrg_time = np.append(avrg_time, middle)

        self.c = avrg_conc
        self.t = avrg_time

        return None

    def fast_find_spectrum(self, f_samps, mag_low, mag_high, raw = False):
        """
        Finds the Lomb-Scargle periodogram of the concentration time series,
        frequency units are in hz. This method uses the Gatspy package, which
        has a fast fft based periodogram calculator.

        f_samps -- int; number of frequency samples to be used in the fitting,
                        they are spaced logarithmically

        mag_low -- float;
        mag_high -- float; the sampled frequency range is
                    [10^mag_low, 10^mag_high]

        raw -- Boolean; if true, self.c_raw is used, otherwise self.c is used

        NB the frequency grid is linear (logarithmic in self.find_spectrum())
        in this method, so many more points should be used.
        """

        assert mag_low < mag_high

        if raw:
            abs_time = self.t_raw
            gas_conc = self.c_raw
        else:
            abs_time = self.t
            gas_conc = self.c

        rel_time = np.array([])
        #the date/datetime objects should be in increasing order, but it should
        #work all the same
        t0 = EPOCH
        for t in abs_time:
            delta_t = t - t0
            rel_time = np.append(rel_time, delta_t.total_seconds())

        fmin = 10 ** mag_low
        fmax = 10 ** mag_high

        df = (fmax - fmin) / f_samps

        #NB uses same x and t in the above section
        model = LombScargleFast().fit(rel_time, gas_conc)
        rel_power = model.score_frequency_grid(fmin, df, f_samps)
        freqs = fmin + df * np.arange(f_samps)
        abs_amplitude = gas_conc.std() * np.sqrt(2*rel_power)

        self.rel_spectrum = np.array([freqs, rel_power])
        self.spectrum = np.array([freqs, abs_amplitude])
        return None

    def plot_ts(self, join_data = False, save = False, raw = False,
                plot_line = False):
        """
        Plots a time series of the concentrations

        join_data -- Boolean; if true, the data are connected with lines.
                     Otherwise only points are plotted with + markers

        save -- Boolean; if true, the figure is saved as a .pdf with
                'self.t_title' as the filename. Otherwise the figure will be
                kept open in the viewer for inspection and must be saved
                manually. There must be a folder called 'Stations' in the
                script directory.

        raw -- if true self.c_raw will be plotted, otherwise self.c will be
               plotted

        plot_line -- Boolean; if true, the 'line' is plotted along with the
                     data. Otherwise it is not plotted.
        """

        if self.line == None:
            self.find_line(raw)
        else:
            pass

        if raw:
            times = self.t_raw
            gas_concs = self.c_raw
        else:
            times = self.t
            gas_concs = self.c

        plt.figure(figsize=(16,10))

        if join_data:
            plt.plot(times, gas_concs)
        else:
            plt.plot(times, gas_concs, ls = '', marker = '+', ms = 1)

        if plot_line:
            plt.plot(times, self.line)
        else:
            pass

        plt.xlabel('Date', fontsize = 18)
        plt.ylabel(self.g_name + ' Concentration / ' + self.g_unit,
                   fontsize = 18)
        plt.title(self.t_title, fontsize = 22)
        if save:
            plt.savefig(self.path + '/' + self.t_name + '.pdf',
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print(self.path + '/' + self.t_name + '.pdf')

        return None

    def plot_ps(self, join_data = False, save = False, days = False,
                  rel = False):
        """
        Plots a spectrum of those produced by 'periodogram'.

        join_data -- Boolean; if true, the data are connected with lines.
                     Otherwise only points are plotted with + markers

        save -- Boolean; if true, the figure is saved as a .pdf with
                'self.f_title' as the filename. Otherwise the figure will be
                kept open in the viewer for inspection and must be saved
                manually. There must be a folder called 'Stations' in the
                script directory.

        days -- Boolean; if true, the amplitudes are plotted against the period
                in days. Otherwise the frequency is plotted in hz

        rel -- Boolean; if true a relative power spectrum is plotted. Otherwise
               an absolute amplitude spectrum is plotted
        """
        if (self.spectrum is None) or (self.rel_spectrum is None):
            print("Call the fast_find_spectrum() method first!")
            return None
        else:
            pass

        if rel:
            spectrum = self.rel_spectrum
            ylab = 'Relative ' + self.g_name + ' \'Power\''
            f_name_new = self.f_name + 'relative'
        else:
            spectrum = self.spectrum
            ylab = self.g_name + ' Amplitude / ' + self.g_unit
            f_name_new = self.f_name

        dep_var = spectrum[1]
        if days:
            indep_var = 1/(spectrum[0]*24*60*60)
            xlab = ['Period', 'days']
        else:
            indep_var = spectrum[0]
            xlab = ['Frequency', 'hz']

        plt.figure(figsize=(16,10))
        if join_data:
            plt.plot(indep_var, dep_var)
        else:
            plt.plot(indep_var, dep_var, ls = '', marker = '+', ms = 1)

        plt.xlabel(xlab[0] + ' / ' + xlab[1], fontsize = 18)
        plt.ylabel(ylab, fontsize = 18)
        plt.title(self.f_title, fontsize = 22)
        if save:
            plt.savefig(self.path + '/' + f_name_new + '.pdf',
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print(self.path + '/' + f_name_new + '.pdf')
        return None

###############################################################################

class WDCGG_TS(AtmosGasTS):
    """
    Defines a child class of AtmosGasTS, with methods used to read in and
    parse the .dat files used by the WDCGG:
    http://ds.data.jma.go.jp/gmd/wdcgg/wdcgg.html
    """
    def __init__(self, station, gas, load_filename):
        """
        Initialises the class based on data found in the .dat files.

        station -- string; the station where the data came from, this will be
                   used in filepaths and plot titles; so make sure the folders
                   exit before run time.

        gas -- string; the gas being examined. An assertion error will be
               thrown if the 'gas' does not match the same gas name in the
               header

        load_filename -- string; the file name/file path + file name where the
                         .dat file is stored. The .dat file must be in the
                         following directory:
                         'Stations/station/gas/load_filename.dat
        """
        self.file_data = []
        self.header = {}

        tot_filepath = 'WDCGG/' + station + '/' + gas + '/'
        tot_filepath += load_filename + '.dat'

        with open(tot_filepath, 'rb') as file:
            for row in file:
                string_row = row.decode()
                if string_row[0] == 'C':
                    #This part depends on colon space (': ') appearing only
                    #to separate 'category: value' in the file header
                    try:
                        key, value = string_row.split(': ')
                        key = key.rstrip('\n')
                        value = value.rstrip('\n')
                        self.header[key] = value
                    except ValueError:
                        print('THIS ROW WAS NOT INCLUDED IN THE HEADER:')
                        print(string_row)
                else:
                    self.file_data.append(string_row.split())

        #converts the file data into more useful data
        times, concs = self.__data()

        #picks out useful parts of the header to complete initialisation
        gas_name = self.header['C18 PARAMETER']
        gas_units = self.header['C21 MEASUREMENT UNIT']



        t_title = gas_name + ' Time Series from ' + station
        f_title = gas_name + ' Periodogram from ' + station
        assert gas_name == gas
        #calls the AtmosGasTS to finish the definition
        super().__init__(gas_name, gas_units, t_title, f_title, times, concs,
                         station, 'WDCGG')
        #The position attributes are assigned to None when super is called
        lat = self.header['C12 LATITUDE']
        long = self.header['C13 LONGITUDE']
        self.set_global_position(lat, long)
        return None

    def __data(self):
        """
        Processes the file data, and returns an array of datetime objects and
        an array concentrations which are mutually ordered. This was written
        only to help the __init__() method, hence the privacy.
        """

        gas_conc = np.array([])
        datetime_array = np.array([])

        for index, row in enumerate(self.file_data):
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
                #line starts at the date in the .dat file, it does not include
                #the header. Commented out; too many lines like these

                #print('Gas concentration is negative; row', index, 'not included:')
                #print(row)
                pass
            else:
                #line starts at the date in the .dat file, it does not include the
                #header
                #print('There is some problem; row', index, 'not included:')
                #print(row)
                pass
        return datetime_array, gas_conc

###############################################################################

class JOIN_TS(AtmosGasTS):
    """
    Defines a child class of AtmosGasTS, with methods used to read in and
    parse the .dat files used by the JOIN website:
    https://join.fz-juelich.de/access/
    """
    def __init__(self, station, gas, load_filename):
        """
        Initialises the class based on data found in the .dat files.

        station -- string; the station where the data came from, this will be
                   used in filepaths and plot titles; so make sure the folders
                   exit before run time.

        gas -- string; the gas being examined. An assertion error will be
               thrown if the 'gas' does not match the same gas name in the
               header

        load_filename -- string; the file name/file path + file name where the
                         .dat file is stored. The .dat file must be in the
                         following directory:
                         'Stations/station/gas/load_filename.dat
        """
        path = 'JOIN/' + station + '/' + gas + '/' + load_filename + '.txt'
        header, times, concs = self.__data(path)
        self.header = header
        #picks out useful parts of the header to complete initialisation
        gas_name = (self.header['#parameter_name']).upper()
        #Check this!!!!!!!!!!!!!!!!!!!
        gas_units = 'ppb'

        t_title = gas_name + ' Time Series from ' + station
        f_title = gas_name + ' Periodogram from ' + station
        assert gas_name == gas

        if self.header['#station_name'] != station:
            warn = 'Given station: ' +  station + ', and file station: '
            warn += self.header['#station_name'] + ', may not be the same.'
            print(warn)
        else:
            pass
        #calls the AtmosGasTS to finish the definition
        super().__init__(gas_name, gas_units, t_title, f_title, times, concs,
                         station, 'JOIN')
        #The position attributes are assigned to None when super is called
        lat = self.header['#station_lat']
        long = self.header['#station_lon']
        #Should be in metres...
        alt = self.header['#station_alt']
        self.set_global_position(lat, long, alt)
        return None

    def __data(self, path):
        """
        Processes the file data, and returns an array of datetime objects and
        an array concentrations which are mutually ordered. This was written
        only to help the __init__() method, hence the privacy.
        """
        header = {}
        datetime_array = np.array([])
        gas_concs = np.array([])


        with open(path, 'r') as tf:
            for index, row in enumerate(tf):
                row = row.replace('\n', '')
                if row[0] == '#':
                    #1 as some ': ' are present after the first one
                    try:
                        row_list = row.split(': ', 1)
                        header[row_list[0]] = row_list[1]
                    except IndexError:
                        print(row, '-- was not included. line:', index+1)
                else:
                    date_and_time, conc = row.split(';')
                    date, time = date_and_time.split(' ')
                    yr, mth, day = date.split('-')
                    hr, mins = time.split(':')

                    yr, mth, day = int(yr), int(mth), int(day)
                    hr, mins = int(hr), int(mins)

                    if float(conc) >= 0:
                        datetime_obj = datetime.datetime(yr, mth, day, hr, mins)
                        datetime_array = np.append(datetime_array, datetime_obj)
                        gas_concs = np.append(gas_concs, float(conc))
                    else:
                        pass
        return header, datetime_array, gas_concs

###############################################################################

class GasCorrelation(object):
    """
    Contains attributes and methods to store the data of two gases, and find
    links between the two data sets.
    """
    def __init__(self, gas_obj0, gas_obj1, interpol_order = 'linear'):
        """
        Stores two gas objects in the attributes. It is assumed that gas_obj0
        has a sampling frequency greater than, or equal too, gas_obj1. If this
        condition is not met, the methods shold still work but there may be
        increased error in interpolation. The gas objects must have some
        overlaping time period for this to work. Also, the interpolated
        concentration is treaded as the independant (i.e. "y" like) variable.

        gas_objX -- instance of a gas class;

        interpol_order -- string (or int); order of the interpolation to be
                          used. See 'AtmosGasTS.set_interpol_data()'.
        """

        self.g0_name = gas_obj0.g_name
        self.g1_name = gas_obj1.g_name
        self.g0_unit = gas_obj0.g_unit
        self.g1_unit = gas_obj1.g_unit

        if gas_obj0.station == gas_obj1.station:
            pass
        else:
            print("The input stations may not be the same")

        self.station = gas_obj0.station
        self.path = gas_obj0.source + '/' + gas_obj0.station + '/'
        self.path += 'Correlations'

        self.t = None
        self.c0 = None
        self.c1 = None
        self.ratio = None
        self.find_interpol_data(gas_obj0, gas_obj1, interpol_order)


        self.t_rel = None
        self.trend_stats = None
        self.line_stats = None
        return None

    def find_interpol_data(self, gas_obj0, gas_obj1, interpol_order):
        """
        This method calls the relevant interpolation methods of the gas objects
        to create concentration attributes in this class.
        """

        #makes sure gas1 covers the same time period as gas 0
        gas_obj1.cut_series(gas_obj0.t)
        #sets data in gas0 to interpolate into the timings of gas1
        gas_obj0.set_interpol_data(gas_obj1.t, order = interpol_order)

        self.t = gas_obj1.t
        self.c0 = gas_obj0.c
        self.c1 = gas_obj1.c
        self.ratio = self.c0/self.c1
        return None

    @staticmethod
    def __t_rel(time_arr):
        """
        Returns an array of relative times (since EPOCH). This method was meant
        for use within other methods, hence privacy.

        time_arr -- array; an array of datetime objects whose relative times
                    will be found.
        """

        rel_time = np.array([])
        t0 = EPOCH
        for t in time_arr:
            delta_t = t - t0
            rel_time = np.append(rel_time, delta_t.total_seconds())
        return rel_time

    def find_t_rel(self):
        """
        Finds the relative times (since the first datetime object).
        """
        self.t_rel = self.__t_rel(self.t)
        return None

    def find_line_statistics(self):
        """
        Finds the statistics of the regresstion line between the concentration
        data. NB the 'line' attribute is an array forming the line of best fit
        of 'c0' as a function of 'c1'.
        """
        reg = linregress(self.c1, self.c0)
        self.line_stats = reg
        self.line = reg[0]*self.c1 + reg[1]
        return None

    def find_trend_statistics(self):
        """
        Finds a line of best fit for the ratio time series as a function of
        time in seconds.
        """

        if self.t_rel == None:
            self.find_t_rel()
        else:
            pass
        trend = linregress(self.t_rel, self.ratio)
        self.trend_stats = trend
        self.trend = trend[0]*self.t_rel + trend[1]
        return None

    def write_reg_data(self, title = 'reg_stats'):
        """
        Writes the regression data into a text file
        """
        print('Regression data saved as:')
        print(self.path + '/' + title + '.txt')
        with open(self.path + '/' + title + '.txt',  'w') as txtfile:
            txtfile.write('These data were measured in ' + self.station + '\n')
            txtfile.write('The units of gradient are, in the first case:\n')
            txtfile.write(self.g0_unit + '/' + self.g1_unit + '\n')
            txtfile.write('In second case:\n')
            txtfile.write(self.g0_unit + '/' + self.g1_unit + '*second')
            txtfile.write('The gases are: ' + self.g0_name + ' and ')
            txtfile.write(self.g1_name + '\n')
            txtfile.write('Regression statistics for the scatter,')
            txtfile.write(' and time series:' + '\n')
            txtfile.write('(gradient, intercept, r value, p value,')
            txtfile.write(' std error in gradient)' + '\n')
            txtfile.write(str(self.line_stats) + '\n')
            txtfile.write(str(self.trend_stats) + '\n')
        return None

    def find_ratio_spectrum(self, f_samps, mag_low, mag_high):
        """
        Finds the Lomb-Scargle periodogram of the concentration time series,
        frequency units are in hz. This method uses the Gatspy package, which
        has a fast fft based periodogram calculator.

        f_samps -- int; number of frequency samples to be used in the fitting,
                        they are spaced logarithmically

        mag_low -- float;
        mag_high -- float; the sampled frequency range is
                    [10^mag_low, 10^mag_high]
        """

        assert mag_low < mag_high
        rel_time = np.array([])
        #the date/datetime objects should be in increasing order, but it should
        #work all the same
        t0 = EPOCH
        for t in self.t:
            delta_t = t - t0
            rel_time = np.append(rel_time, delta_t.total_seconds())

        fmin = 10 ** mag_low
        fmax = 10 ** mag_high

        df = (fmax - fmin) / f_samps

        #NB uses same x and t in the above section
        model = LombScargleFast().fit(rel_time, self.ratio)
        rel_power = model.score_frequency_grid(fmin, df, f_samps)
        freqs = fmin + df * np.arange(f_samps)
        abs_amplitude = self.ratio.std() * np.sqrt(2*rel_power)

        self.rel_spectrum = np.array([freqs, rel_power])
        self.spectrum = np.array([freqs, abs_amplitude])
        return None

    def plot_ratio_spectrum(self,  title = '', fig_name = '',
                            join_data = False, save = False, days = False,
                            rel = False):
        """
        Plots a spectrum of those produced by 'periodogram'.

        title -- string; title that will appear on the plot.

        fig_name -- string; will be the filename used if 'save' is true.

        join_data -- Boolean; if true, the data are connected with lines.
                     Otherwise only points are plotted with + markers

        save -- Boolean; if true, the figure is saved as a .pdf with
                'self.f_title' as the filename. Otherwise the figure will be
                kept open in the viewer for inspection and must be saved
                manually. There must be a folder called 'Figures' in the script
                directory.

        days -- Boolean; if true, the amplitudes are plotted against the period
                in days. Otherwise the frequency is plotted in hz

        rel -- Boolean; if true a relative power spectrum is plotted. Otherwise
               an absolute amplitude spectrum is plotted
        """
        if (title != '') and (fig_name == ''):
            fig_name = title
        else:
            pass

        if (self.spectrum == None) or (self.rel_spectrum == None):
            print("Call the get_ratio_spectrum() method first!")
            return None
        else:
            pass

        if rel:
            spectrum = self.rel_spectrum
            ylab = 'Ratio Relative Power'
        else:
            spectrum = self.spectrum
            ylab = 'Ratio Amplitude'
            if self.g0_unit == self.g1_unit:
                pass
            else:
                ylab += self.g0_unit + '/' + self.g1_unit

        dep_var = spectrum[1]
        if days:
            indep_var = 1/(spectrum[0]*24*60*60)
            xlab = ['Period', 'days']
        else:
            indep_var = spectrum[0]
            xlab = ['Frequency', 'hz']

        plt.figure(figsize=(16,10))
        if join_data:
            plt.plot(indep_var, dep_var)
        else:
            plt.plot(indep_var, dep_var, ls = '', marker = '+', ms = 1)

        plt.xlabel(xlab[0] + ' / ' + xlab[1], fontsize = 18)
        plt.ylabel(ylab, fontsize = 18)
        plt.title(title, fontsize = 22)
        if save:
            plt.savefig(self.path + '/' + fig_name + '.pdf',
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print(self.path + '/' + fig_name + '.pdf')
        return None

    def plot_ratio_ts(self, title = '', fig_name = '', join_data = False,
                      save = False, trend = True):
        """
        Plots the time series of the ratio of c0 to c1.

        title -- string; title that will appear on the plot.

        fig_name -- string; will be the filename used if 'save' is true.

        join_data -- Boolean; if true, the data are connected with lines.
                     Otherwise only points are plotted with + markers.

        join_data -- Boolean; if true, the data are connected with lines.
                 Otherwise only points are plotted with + markers.

        save -- Boolean; if true, the figure is saved as a .pdf with
                'fig_name' as the filename. Otherwise the figure will be
                kept open in the viewer for inspection and must be saved
                manually. There must be a folder called 'Figures/GasCorrs'
                in the script directory.

        trend -- Boolean; if true, the trendline of the ratio is plotted also.
                 Otherwise the trendline is not plotted
        """
        if (title != '') and (fig_name == ''):
            fig_name = title
        else:
            pass

        ylab = 'Ratio of ' +  self.g0_name + ' to ' + self.g1_name
        if self.g0_unit == self.g1_unit:
            pass
        else:
            ylab += ' / ' + self.g0_unit + '/' + self.g1_unit


        plt.figure(figsize=(16,10))
        if join_data:
            plt.plot(self.t, self.ratio)
        else:
            plt.plot(self.t, self.ratio, ls = '', marker = '+', ms = 2)

        if trend:
            plt.plot(self.t, self.trend)
        else:
            pass

        plt.xlabel('Date', fontsize = 18)
        plt.ylabel(ylab, fontsize = 18)
        plt.title(title, fontsize = 22)

        if save:
            plt.savefig(self.path + '/' + fig_name + '.pdf',
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print(self.path + '/' + fig_name + '.pdf')
        else:
            pass
        return None

    def plot_scatter(self, title = '', fig_name = '', save = False,
                     line = True):
        """
        Plots a scatter chart of the interpolated concentration data of 'c0'
        against the uninterpolated concentraion data 'c1'.

        title -- string; title that will appear on the plot.

        fig_name -- string; will be the filename used if 'save' is true.

        join_data -- Boolean; if true, the data are connected with lines.
                     Otherwise only points are plotted with + markers.

        save -- Boolean; if true, the figure is saved as a .pdf with
                'fig_name' as the filename. Otherwise the figure will be
                kept open in the viewer for inspection and must be saved
                manually. There must be a folder called 'Figures/GasCorrs'
                in the script directory.
        """
        if (title != '') and (fig_name == ''):
            fig_name = title
        else:
            pass

        plt.figure(figsize=(16,10))
        plt.plot(self.c1, self.c0, ls = '', marker = '+', ms = 4)
        if line:
            plt.plot(self.c1, self.line)
        else:
            pass

        plt.xlabel(self.g1_name + ' / ' + self.g1_unit, fontsize = 18)
        plt.ylabel(self.g0_name + ' / ' + self.g0_unit, fontsize = 18)
        plt.title(title, fontsize = 22)

        if save:
            plt.savefig(self.path + '/' + fig_name + '.pdf',
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print(self.path + '/' + fig_name + '.pdf')
        else:
            pass
        return None

###############################################################################

def ts_multi_plot(gas_obj_arr, title = '', fig_name = '', join_data = False,
                  save = False, raw = False):
    """
    Takes an array of gas objects, and plots the time series on one plot. We
    will assume that all the gas concentrations are in ppb in all passed
    classes.

    gas_obj_arr -- array like; iterable of multiple gas (eg AtmosGasTS)
                   objects. They must be fully initialised and processed.

    title -- string; title that will appear on the plot.

    fig_name -- string; will be the filename used if 'save' is true.

    join_data -- Boolean; if true, the data are connected with lines.
                 Otherwise only points are plotted with + markers.

    save -- Boolean; if true, the figure is saved as a .pdf with
            'fig_name' as the filename. Otherwise the figure will be
            kept open in the viewer for inspection and must be saved
            manually. There must be a folder called 'Figures/MultiPlots'
            in the script directory.

    raw -- Boolean; if true the raw data in each class is plotted. Otherwise
           the processed data are plotted.
    """

    plt.figure(figsize=(16,10))

    for gas_obj in gas_obj_arr:

        if raw:
            times = gas_obj.t_raw
            gas_concs = gas_obj.c_raw
        else:
            times = gas_obj.t
            gas_concs = gas_obj.c

        lab = gas_obj.station + ' ' + gas_obj.g_name
        if join_data:
            plt.plot(times, gas_concs, label = lab)
        else:
            plt.plot(times, gas_concs, ls = '', marker = '+', ms = 1,
                     label = lab)

    plt.xlabel('Date', fontsize = 18)
    plt.ylabel('Concentration / ppb', fontsize = 18)

    plt.title(title, fontsize = 22)
    plt.legend(loc = 'best')

    if save:
        plt.savefig('Figures/MultiPlots/' + fig_name + '.pdf', format = 'pdf')
        plt.close()
        print('Figure saved as:')
        print('Figures/MultiPlots/' + fig_name + '.pdf')

    return None


def ps_multi_plot(gas_obj_arr, title = '', fig_name = '', join_data = False,
                  save = False, days = False):
        """
        Plots the spectra of the gases in gas_obj_arr. Assumed units are ppb.

        gas_obj_arr -- array like; iterable of multiple gas (eg AtmosGasTS)
                   objects. They must be fully initialised and processed.

        join_data -- Boolean; if true, the data are connected with lines.
                     Otherwise only points are plotted with + markers

        title -- string; title that will appear on the plot.

        fig_name -- string; will be the filename used if 'save' is true.

        save -- Boolean; if true, the figure is saved as a .pdf with
                'self.f_title' as the filename. Otherwise the figure will be
                kept open in the viewer for inspection and must be saved
                manually. There must be a folder called 'Figures/MultiPlots'
                in the script directory.

        days -- Boolean; if true, the amplitudes are plotted against the period
                in days. Otherwise the frequency is plotted in hz
        """

        plt.figure(figsize=(16,10))

        for gas_obj in gas_obj_arr:
            if (gas_obj.spectrum == None) or (gas_obj.rel_spectrum == None):
                print("Call the fast_find_spectrum() method first!")
                return None
            else:
                pass

            spectrum = gas_obj.spectrum
            dep_var = spectrum[1]
            if days:
                indep_var = 1/(spectrum[0]*24*60*60)
                xlab = ['Period', 'days']
            else:
                indep_var = spectrum[0]
                xlab = ['Frequency', 'hz']

            lab = gas_obj.station + ' ' + gas_obj.g_name
            if join_data:
                plt.plot(indep_var, dep_var, label = lab)
            else:
                plt.plot(indep_var, dep_var, ls = '', marker = '+', ms = 1,
                         label = lab)

        ylab = 'Amplitude / ppb'
        plt.xlabel(xlab[0] + ' / ' + xlab[1], fontsize = 18)
        plt.ylabel(ylab, fontsize = 18)

        plt.title(title, fontsize = 22)
        plt.legend(loc = 'best')
        if save:
            plt.savefig('Figures/MultiPlots' + fig_name + '.pdf',
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print('Figures/MultiPlots' + fig_name + '.pdf')
        return None
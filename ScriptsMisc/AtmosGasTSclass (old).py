# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:12:56 2016

@author: tmed2
"""


import datetime
import numpy as np
import matplotlib.pyplot as plt
#may remove the scipy based periodogram finder
from scipy import signal
from scipy.stats import linregress
from scipy.interpolate import interp1d
from gatspy.periodic import LombScargleFast


#This is used to define the start of time, so that relative times will be
#consistent
EPOCH = datetime.datetime(1975, 1, 1)


class AtmosGasTS(object):
    """
    Class containing attributes and methods relating to the time series
    of atmospheric gases.
    """
    def __init__(self, gas_name, gas_units, t_title, f_title, times, concs,
                 subpath = ''):
        """
        Initialises the AtmosGasClass.
        
        gas_name -- string;
        
        gas_units -- string;
        
        t_title -- string; name to be used on time series.
        
        f_title -- string; name to be used on power spectra.
        
        times -- array; an array of date/datetime objects ordered with concs.
        
        concs -- array; an array of the gas concentrations corresponding to
                 the 'times' array.
        
        subpath -- string; file path of saved figures: Figures/subpath/title.pdf
                 Default is an empty string (which saves in Figures/title.pdf).
                 Intened to be one folder deep; the '/' should not be included 
                 unless there is more than one folder in the middle. The path 
                 must exist before the scripts are run.
        """
               
        
        self.g_name = gas_name
        self.g_unit = gas_units
        self.t_title = t_title
        self.f_title = f_title
        
        self.path = subpath
        #these will be used to name the files, they are the same as x_title
        #but with newlines ('\n') replaced with spaces (' ')
        self.t_name = t_title.replace('\n', ' ')
        self.f_name = f_title.replace('\n', ' ')
        
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
        if self.t_rel == None:
            self.find_t_rel(raw)
        else:
            pass
        
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
        
        if self.line_stats == None:
            self.find_line_stats(raw)
        else:
            pass
        
        gradient = self.line_stats[0]
        intercept = self.line_stats[1]
        
        self.line = intercept + self.t_rel * gradient
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
        self.find_t_rel(raw = raw)
        if raw:
            times = self.t_rel
            concs = self.c_raw
        else:
            times = self.t_rel
            concs = self.c
        
        self.interpol = interp1d(times, concs, kind = order)
        return None
    
    def set_interpol_data(self, time_arr, raw = False):
        """
        Uses the 'interpol' method-attribute to find the value of the function
        at the points in time_arr. This method overwrites the the processed
        data attributes.
        
        time_arr -- array; array of times (including datetime objects) that is
                    used to calculate the interpolated data.
        
        raw -- Boolean; if True the raw times are used.
        """
            
        if raw:
            t = self.t_raw
        else:
            t = self.t
        low = np.min(t)
        high = np.max(t)
        
        if self.interpol == None:
            self.find_interpol(raw = raw)
        else:
            pass
        
        #this ensure that 'time_arr' is within the range of interpolation
        time_arr = time_arr[(time_arr > low) & (time_arr < high)]
        self.t = time_arr
        
        time_arr = self.__t_rel(time_arr)
        self.c = self.interpol(time_arr)
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
        average_gas_conc = np.mean(self.c)
        sigma_gas_conc = np.std(self.c)
        upper_cut = average_gas_conc + sigma_cut*sigma_gas_conc
        lower_cut = average_gas_conc - sigma_cut*sigma_gas_conc
        
        new_gas_conc = np.array([])
        new_time = np.array([])
        
        #vectorise this at some point...
        for index, conc in enumerate(self.c):
            if (conc <= upper_cut) and (conc >= lower_cut):
                new_gas_conc = np.append(new_gas_conc, conc)
                new_time = np.append(new_time, self.t_raw[index])
                
            elif (conc > upper_cut) or (conc < lower_cut):
                #this concerns the input time series, index starts at 0
                #print('Index ', index, 
                      #'removed from c_raw; concentration out ofrange')
                pass
            else:
                #this concerns the input time series, index starts at 0
                #print('Index ', index, 
                      #'removed from c_raw; there is some problem')
                pass
        
        self.t = new_time
        self.c = new_gas_conc
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
        
    #NB, the default values may be different to functions in my other scripts!
    def ts_plot(self, join_data = False, save = False, raw = False, 
                plot_line = False):
        """ 
        Plots a time series of the concentrations
        
        join_data -- Boolean; if true, the data are connected with lines. 
                     Otherwise only points are plotted with + markers
                     
        save -- Boolean; if true, the figure is saved as a .pdf with 
                'self.t_title' as the filename. Otherwise the figure will be 
                kept open in the viewer for inspection and must be saved 
                manually. There must be a folder called 'Figures' in the script
                directory.
                
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
            plt.savefig('Figures/' + self.path + '/' + self.t_name + '.pdf', 
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print('Figures/' + self.path + '/' + self.t_name + '.pdf')
            
        return None
    
    #Use fast_find_spectrum instead
    def find_spectrum(self, f_samps, mag_low, mag_high, raw = False):
        """
        Finds the Lomb-Scargle periodogram of the concentration time series,
        frequency units are in hz.
        
        f_samps -- int; number of frequency samples to be used in the fitting,
                        they are spaced logarithmically
                        
        mag_low -- float;
        mag_high -- float; the sampled frequency range is 
                    [10^mag_low, 10^mag_high]
        
        raw -- Boolean; if true, self.c_raw is used, otherwise self.c is used
        """
        
        assert mag_low < mag_high
        
        if raw:            
            abs_time = self.t_raw
            gas_conc = self.c_raw
        else:
            abs_time = self.t
            gas_conc = self.c
        
        rel_time = np.array([])
        angfreqs = 2*np.pi*np.logspace(mag_low, mag_high, f_samps)
        #the date/datetime objects should be in increasing order, but it should
        #work all the same
        t0 = abs_time[0]
        for t in abs_time:
            delta_t = t - t0
            rel_time = np.append(rel_time, delta_t.total_seconds())
            
        rel_gas_conc = gas_conc - np.mean(gas_conc)
        
        pgram = signal.lombscargle(rel_time, rel_gas_conc, angfreqs)
        N = rel_time.shape[0]
        norm_pgram = np.sqrt(4*(pgram/N))
        
        norm_pow_spectrum = np.array([angfreqs/(2*np.pi), norm_pgram])
        self.spectrum = norm_pow_spectrum
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
        
    def ps_plot(self, join_data = False, save = False, days = False, 
                  rel = False):
        """
        Plots a spectrum of those produced by 'periodogram'. 
        
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
        if (self.spectrum == None) or (self.rel_spectrum == None):
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
            plt.savefig('Figures/' + self.path + '/' + f_name_new + '.pdf', 
                        format = 'pdf')
            plt.close()
            print('Figure saved as:')
            print('Figures/' + self.path + '/' + f_name_new + '.pdf')
        return None
        
        
class WDCGG_TS(AtmosGasTS):
    """
    Defines a child class of AtmosGasTS, with methods used to read in and
    parse the .dat files used by the WDCGG:
    http://ds.data.jma.go.jp/gmd/wdcgg/wdcgg.html
    """
    def __init__(self, country, load_filename):
        """
        Initialises the class based on data found in the .dat files.
        
        country -- string; the nation where the data came from, this will be
                   used in filepaths and plot titles; so make sure the folders 
                   exit before run time.
        
        load_filename -- string; the file name/file path + file name where the
                         .dat file is stored. The .dat file must be in the 
                         following directory:
                         'Data/WDCGG-SurfaceData/country/.../load_filename.dat
                         or in a sub direcory of this which can (and should) be
                         included in 'load_filename'.
        """
        self.file_data = []
        self.header = {}
        
        tot_filepath = 'Data/WDCGG-SurfaceData/' + country + '/'
        tot_filepath +=load_filename + '.dat'
        
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
        self.station = self.header['C07 STATION NAME']
        
        gas_name = self.header['C18 PARAMETER']
        gas_units = self.header['C21 MEASUREMENT UNIT']
        t_title = gas_name + ' Time Series from ' + self.station
        f_title = gas_name + ' Periodogram from ' + self.station
        spath = 'WDCGG/' + country +'/' + gas_name
        
        #calls the AtmosGasTS to finish the definition
        super().__init__(gas_name, gas_units, t_title, f_title, times, concs,
                         subpath = spath)
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
            'self.f_title' as the filename. Otherwise the figure will be 
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
            
        
        if join_data:
            plt.plot(times, gas_concs, label = gas_obj.g_name)
        else:
            plt.plot(times, gas_concs, ls = '', marker = '+', ms = 1, 
                     label = gas_obj.g_name)
            
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
            
          
            if join_data:
                plt.plot(indep_var, dep_var, label = gas_obj.g_name)
            else:
                plt.plot(indep_var, dep_var, ls = '', marker = '+', ms = 1,
                         label = gas_obj.g_name)
                
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

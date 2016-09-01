# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:25:14 2016

@author: tmed2
"""


import datetime
import numpy as np
import matplotlib.pyplot as plt
import AtmosGasTools.AtmosGasTS as ag


def comparison_wrapper(station, g0name, g1name, g0_file, g1_file,
                       all_save = False, all_join = False, dev_check = 0):
    """
    Acts as a wrapper for the comparison class and its methods. g0 should have
    an equal or greater sampling rate than g1, as g0 is being interpolated
    using the timings of g1.

    station -- string; the physcial site of the measurements, will be used
               in the file directories

    gXname -- string; name of the Xth gas, will be used in file directories

    gX_file -- string; name of the file containing the relevant data of the Xth
               gas, don't include the file extension (.dat).

    all_save -- Boolean; if true, all the figures will be saved as will the
                regression data.

    all_join -- Boolean; if true, all the figureswill have consecutive points
                connected

    dev_check -- float; if > 0, points that are (two sided) dev_check standard
                 deviations away from the mean will be removed
    """
    g0 = ag.WDCGG_TS(station, g0name, g0_file)
    g1 = ag.WDCGG_TS(station, g1name, g1_file)

    if dev_check > 0:
        g0.deviation_check(dev_check)
        g1.deviation_check(dev_check)
    else:
        pass
    
    #the regression data are found, implicitly, in the 'plot_ts()' method
    g0.plot_ts(save = all_save, join_data = all_join, plot_line = True)
    g1.plot_ts(save = all_save, join_data = all_join, plot_line = True)
    
    if all_save:
         g0.write_reg_data()
         g1.write_reg_data()
    else:
        pass

    corr_obj = ag.GasCorrelation(g0, g1, interpol_order = 'linear')
    corr_obj.get_trend_statistics()
    corr_obj.get_line_statistics()

    #NB this may have to changed depending on the data
    corr_obj.get_ratio_spectrum(250000, -8.0, -3.85)

    title_ps = 'Spectrum of the Ratio of ' + g0name + ' to ' + g1name
    title_ps += ' from ' + station
    corr_obj.plot_ratio_spectrum(title = title_ps, save = all_save,
                             join_data = all_join, days = True)

    title_ts = 'Ratio of ' + g0name + ' to ' + g1name + ' from ' + station
    corr_obj.plot_ratio_ts(title = title_ts, save = all_save, join_data = all_join)

    title_sc = g0name + ' as a Function of '+ g1name + ' from ' + station
    corr_obj.plot_scatter(title = title_sc, save = all_save)
    
    if all_save:
        corr_obj.write_reg_data()
    else:
        pass
    
    return None

#uncomment if needed
#comparison_wrapper('Heimaey', 'O3', 'CO', 'Iceland O3 2003-2010',
#                   'HeimaeyCO', all_save = True, all_join = True)
#
#comparison_wrapper('Izana', 'O3', 'CO', 'Izana O3 2008-2014',
#                   'Izana CO 2008-2014', all_save = True, all_join = True)
#
#comparison_wrapper('Mace Head', 'O3', 'CO', 'MaceHeadO3', 'MaceHeadCO',
#                   all_save = True, all_join = True)
#
#comparison_wrapper('Pic du Midi', 'O3', 'CO', 'Pic du Midi O3 2008-2012',
#                   'Pic du Midi CO 2008-2012', all_save = True,
#                   all_join = True, dev_check = 4)
#
#comparison_wrapper('Tudor Hill', 'O3', 'CO', 'TudorHillO3 1989-1998',
#                   'TudorHillCOMonth', all_save = True, all_join = True)

comparison_wrapper('Ragged Point', 'O3', 'CO', 'BarbadosO3 2006-2015',
                   'Ragged Point CO', all_save = True, all_join = True)











#this was the base of the 'comparison_wrapper()' function
"""
###############################################################################
#This is  O3-CO analysis of Mace Head, Ireland
O3 = ag.WDCGG_TS('Izana', 'O3', 'Izana O3 2008-2014')
CO = ag.WDCGG_TS('Izana', 'CO', 'Izana CO 2008-2014')

all_save = True

CO.plot_ts(save = all_save, raw = True)
O3.plot_ts(save = all_save, raw = True)

both = ag.GasCorrelation(O3, CO, interpol_order = 'linear')
both.get_trend_statistics()
both.get_line_statistics()
both.get_ratio_spectrum(250000, -8.0, -3.85)


both.plot_ratio_spectrum(title = 'Spectrum of the Ratio of Ozone to Carbon' + \
                         ' Monoxide from Izana', save = all_save,
                         join_data = True, days = True)

both.plot_ratio_ts(title = 'Ratio of Ozone to Carbon Monoxide from Izana',
                   save = all_save, join_data = True)

both.plot_scatter(title = 'Ozone as a function of Carbon Monoxide from' + \
                  ' Izana', save = all_save)
both.write_reg_data()
###############################################################################
"""
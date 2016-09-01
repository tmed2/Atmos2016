# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 13:47:29 2016

@author: tmed2

This module provides functions to join hourly data files (from the WDCGG:
http://ds.data.jma.go.jp/gmd/wdcgg/) into one large file that can be processed
as normal.

NB:
#File paths must be relative to the calling script (not this, the defining
 script)

Available functions:
-dat_reader: reads in an individual hourly data file
-bin_writer: writes input data into a new file, in append mode
-data_tree_join: combines all the hourly data files into one file
"""


def dat_reader(fpath, fname):
    """
    Reads in an individual file
    """

    header = []
    data = []
    with open(fpath + fname + '.dat', 'rb') as file:
        for row in file:
            string_row = row.decode('iso-8859-1')
            if string_row[0] == 'C':
                header.append(string_row)
            else:
                data.append(string_row)

    return [header, data]


def bin_writer(fpath, fname, data):
    """
    Writes 'data' onto the end of a file (or creates it if it does not exist).
    """
    path = fpath + fname + '.dat'
    with open(path, 'ab') as file:
        for row in data:
            file.write(row.encode('utf-8'))
    return None



def data_tree_join(fpath, conv, ystart, yend, joined_fname):
    """
    Looks at the yearly folders in 'fpath' to pull the data together, writing a
    new file of joined data in the same location.

    fpath -- string; the directory of the yXXXX year folders, there must be
                     a '/' on the end

    conv -- string; the convention used to name the files. Assumed to all be
                    of the form "convXXXX" where XXXX is the year (don't
                    include the year, or file ending). e.g. the hourly data
                    for O3 in 1989 is "rpb413n00.noaa.as.cn.o3.nl.hr1989.dat",
                    so the 'conv' is "rpb413n00.noaa.as.cn.o3.nl.hr"

    ystart -- int; year to start with (4 digit format)

    yend -- int; year to end with (4 digit format)

    joined_fname -- string; the name of the file containg the joined_data
    """

    head_write = True
    for year in range(ystart, yend + 1):
        file_loc = conv + str(year)
        fpath_loc = fpath + 'y' + str(year) + '/'

        try:
            header, data = dat_reader(fpath_loc, file_loc)
            if head_write:
                #The header of the first valid file is used, the user must be
                #sure that the value does not change across headers
                bin_writer(fpath, joined_fname, header)
                head_write = False
            else:
                pass
            bin_writer(fpath, joined_fname, data)
        except FileNotFoundError:
            print('No data found for', year)
    return None

"""
#Example usage
data_tree_join('SplitData/CapeVerdeO3Hourly/', 'cvo116n00.uyrk.as.cn.o3.nl.hr',
               2006, 2015, 'Cape Verde O3 2006-2015')
"""

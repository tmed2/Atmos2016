###############################################################################
The 'atmosgastools' package contains various tools used to process, analyse,
and display time series data of atmospheric gases. It consists of the following
modules:
---atmosgasts
---wdcggtools
---nasatools
More detailed information can be found in the headers of each of the module's
files, and their code.

THE FOLLOWING IS ESSENTIAL:
File paths must be relative to the calling script, as opposed to the modules
which are defining scripts.

In order to use the the AtmosGasTools package the following directory structure
must be present. <category> is a placeholder for a specific case of that
category; # precedes comments; ... = etc.

\AtmosGasTools
	\__init__.py
	...
	
\<source1> #e.g. JOIN, WCDGG, NASA...
	\<station> #where the data was collected
		\<gas>
			\<data>.<file_type> #e.g. O3_hourly.dat
		\Correlations

\<source2>
	...

\SplitData
	\<station>
		\y<year1>
		\y<year1>
		...

\Figures
	\MultiPlots #multiplot function figures are saved here

\<main_script1>.py #calls AtmosGasTools and it's modules
\<main_script2>.py
...
\<other_script1>.py #does not call AtmosGasTools
...
###############################################################################
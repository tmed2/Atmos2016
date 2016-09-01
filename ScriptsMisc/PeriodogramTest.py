# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 09:41:02 2016

@author: tmed2

Example use of the Lomb-Scargle periodogram (LSP):
http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lombscargle.html

Examples for the gatspy LSP:
https://jakevdp.github.io/blog/2015/06/13/lomb-scargle-in-python/
https://github.com/astroML/gatspy/
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
###############################################################################
from scipy import signal

start0 = datetime.datetime.now()

ang = 2*np.pi
t = np.linspace(0, 4, 20000)
x = np.sin(ang*t) + 4*np.cos(5*ang*t) + 5

#these are the sample angular frequencies that the lombscargle function
#will attempt to fit
ang_f = ang*np.linspace(0.001, 10, 1000)

"""
NB the input frequency is angular, and the mean must be subtracted from the
input x values in general (if it is not already zero). otherwise you'll get a
junk spectrum, particularly at the low frequencies
"""
x_tilde = signal.lombscargle(t, x, ang_f)
N = x.shape[0]
#lombscargle returns intensity values. We, might, want amplitude values
x_tilde_norm = 2*np.sqrt((x_tilde/N))

#the time domain
plt.figure(0)
plt.plot(t, x)

#the frequency domain, component frequencies are clear (in hz)
plt.figure(1)
plt.plot(ang_f/ang, x_tilde_norm)

end0 = datetime.datetime.now()
print("SciPy took:", end0 - start0)
###############################################################################
from gatspy.periodic import LombScargleFast

start1 = datetime.datetime.now()

N = 1000
fmin = 0.001
fmax = 10

df = (fmax - fmin) / N

#NB uses same x and t in the above section
model = LombScargleFast().fit(t, x)
power = model.score_frequency_grid(fmin, df, N)
freqs = fmin + df * np.arange(N)

# plot the results
plt.figure()
plt.plot(freqs, ((x.std())*np.sqrt(2*power)),  ls = '', marker = '+')
       
end1 = datetime.datetime.now()
print("Gatspy took:", end1 - start1)
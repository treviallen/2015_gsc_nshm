# -*- coding: utf-8 -*-
"""
Created on Tue Sep 06 10:31:49 2016

Calculates PGV for a given mag and distance for a specified lon/lat

@author: tallen
"""

from sys import argv
from numpy import loadtxt, savetxt, hstack, zeros_like, sqrt, exp, nan
from oq_tools import distance
#from calc_oq_gmpes import atkinson_adams2013_gsim, ena_gsims
from multiprocessing import Pool
import time

t = time.time()
print t

# get lat/lon
loclo = float(argv[1])
locla = float(argv[2])

# set ses file
sesfile = 'SE_HY_event_set.1000000.csv'

# first parse ses file
print 'parsing data...'
data = loadtxt(sesfile, delimiter=',')
#lola = hstack((evlo, evla))
evmw = data[:,0]
evpr = data[:,1]
evlo = data[:,2]
evla = data[:,3]
evdp = data[:,4]

# set vs30 for imt calculation
vs30 = 450
hdf5file = '/home/tallen/oq-engine/2015_gsc_nshm/gm_tables/ENA_med_clC.hdf5'

# make data dict
print 'making data dict...'
datadict = []
for i in range(0, len(evmw)):
    datadict.append({'evmw':evmw[i], 'evlo':evlo[i], 'evla':evla[i], 'evdp':evdp[i]})
    

##############################################################################
# function for multiporcessing
from calc_oq_gmpes import atkinson_adams2013_gsim, ena_gsims
from scipy.constants import g
    
def get_ims(tdict):
    # extract from dict
    evla = tdict['evla']
    evlo = tdict['evlo']
    evdp = tdict['evdp']
    evmw = tdict['evmw']
    
    # calulate distance in km   
    edist = distance(locla, loclo, evla, evlo)[0]
    hdist = sqrt(edist**2 + evdp**2)
    rjb = edist
    rrup = hdist
    
    # get IMs
    if edist < 790.:
        aa13imt = atkinson_adams2013_gsim(evmw, evdp, evdp, rrup, rjb, hdist, hdf5file)
        aa13 = exp(aa13imt['pgv'][0][0]) / g # wierd relic for GSCFRISK
        
        # get AB06 & YA15
        dip = 30.
        rke = 90.
        AB06imt, YA15imt = ena_gsims(evmw, evdp, evdp, dip, rke, rrup, rjb, vs30)
        ab06 = exp(AB06imt['pgv'][0][0])
        ya15 = exp(YA15imt['pgv'][0][0])
        
        #evdi = edist
        tdict['ab06'] = ab06
        tdict['aa13'] = aa13
        tdict['ya15'] = ya15

    else:
        tdict['ab06'] = nan
        tdict['aa13'] = nan
        tdict['ya15'] = nan
    tdict['evdi'] = edist

    return tdict

##############################################################################

#print get_ims(tdict)

#print get_ims(datadict[0])
print 'getting ground motions...'
print '    starting multiprocessing...'
p = Pool()
datadict = p.map(get_ims, datadict)

# prep arrays to be filled
evdi = zeros_like(evmw)
ab06 = zeros_like(evmw)
aa13 = zeros_like(evmw)
ya15 = zeros_like(evmw)

# fill new arrays
for i, dd in enumerate(datadict):
    evdi[i] = dd['evdi']
    ab06[i] = dd['ab06']
    aa13[i] = dd['aa13']
    ya15[i] = dd['ya15']

# refmt data
evmw = data[:,0].reshape(len(data[:,0]), 1)
evpr = data[:,1].reshape(len(data[:,0]), 1)
evlo = data[:,2].reshape(len(data[:,0]), 1)
evla = data[:,3].reshape(len(data[:,0]), 1)
evdp = data[:,4].reshape(len(data[:,0]), 1)

evdi = evdi.reshape(len(evdi), 1)
ab06 = ab06.reshape(len(evdi), 1)
aa13 = aa13.reshape(len(evdi), 1)
ya15 = ya15.reshape(len(evdi), 1)

        
# concat datat
new_data = hstack((data, evdi, ab06, aa13, ya15))

# write to file
outfile = sesfile[0:-3]+'_'.join((str(loclo), str(locla), 'pgv.csv'))
header = 'MAG,OCC_RATE,LON,LAT,DEPTH,EDIST,AB06,AA13,YA15'
savetxt(outfile, new_data, fmt='%.5e', delimiter=',', header=header)

print time.time() - t


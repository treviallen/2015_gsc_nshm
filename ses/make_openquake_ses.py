# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:52:14 2015

Code generates stochastic event set based on OQ source model

Usage:
    python make_openquake_ses.py <sourcexml>
    
    argv sourcexml: path to source model

@author: tallen
"""
from openquake.commonlib.sourceconverter import SourceConverter
from openquake.commonlib.source import SourceModelParser
from openquake.hazardlib.calc import stochastic_event_set
from datetime import datetime
from sys import argv

time = datetime.now()
print 'Start Time:', time.strftime('%Y-%m-%d %H:%M:%S')

itime = 1000000.0
#itime = 200.0
discarea = 10.0
converter = SourceConverter(itime, # Investigation time = Length of Synthetic Catalogue
                            5.0, # Simple fault rupture mesh spacing
                            5.0, # complex fault rupture mesh spacing
                            0.1, # Width of MFD bin
                            discarea # Area source discretisation
                            )

sourcexml = argv[1]
src_model = SourceModelParser(converter)
groups = src_model.parse_groups(sourcexml)


# put sources back into a single list
sources = []
for group in groups:
    for src in group:
        sources.append(src)

ses = list(stochastic_event_set(sources))
sesmag = []
sesdep = []
seslat = []
seslon = []

# export mag, hypocenter, rake and surface for each rupture
print '\nExporting Event Set...\n'
f = open("event_set.csv", "w")
for rup in ses:
    print >> f, "{:s}".format(",".join([str(val) for val in [rup.mag, \
                              rup.occurrence_rate, rup.hypocenter.longitude, \
                              rup.hypocenter.latitude, rup.hypocenter.depth]]))
    '''
    sesmag.append(rup.mag)
    sesdep.append(rup.hypocenter.depth)
    seslon.append(rup.hypocenter.longitude)
    seslat.append(rup.hypocenter.latitude)
    
    # get ruptuures
    rup.surface.get_bounding_box()
    rup.surface.bottom_left.depth
    rup.surface.top_left.depth
    '''

        
f.close()
time = datetime.now()
print 'Stop Time:', time.strftime('%Y-%m-%d %H:%M:%S'), '\n'
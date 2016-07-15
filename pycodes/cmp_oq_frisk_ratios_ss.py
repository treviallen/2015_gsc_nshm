from __future__ import unicode_literals
from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
#from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
from numpy import array, diff, exp, log, interp, loadtxt, vstack, mean
from oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path, sep
import warnings, sys
reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

'''
Calculates ratios between OQ and GSCFRISK for sigle sources only

file paths must be manually specified
'''
###############################################################################
# settings
###############################################################################

period = '0.2'
job = 'chv_nbcc'
job_num = '601'

tmpper = period.replace('.','')

# set frisk files
friskfolder = 'CHVcomparison'
friskfile   = 'OQ_CHVregion_0.2_20s_000thperc.sol'


'''
job = sys.argv[1]
job_num = sys.argv[2]
period = sys.argv[3]
'''

'''
   410 | successful | 2016-05-20 21:20:26 | SWCan 2015 NBCC Hazard - ALL BRANCHES - swcan # not complex faults
   416 | successful | 2016-05-21 00:13:58 | SECan 2015 NBCC Hazard - ALL BRANCHES - secan
   
'''   
   


###############################################################################
# read frisk data
###############################################################################

friskpath = path.join('','data','gscfrisk_comparisons',friskfolder,friskfile)
friskhaz = loadtxt(friskpath, delimiter=',', skiprows=4, usecols = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
friskprobs = array([0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001])

# get place names
places = []
lines = open(friskpath).readlines()[4:]
for line in lines:
    places.append(', '.join(line.strip().split(',')[-2:]))


###############################################################################
# read OQ data
###############################################################################
hazcurvefile = path.join('','jobs_hazard',job,'out','hazard_curve-mean_'+job_num+'-SA('+period+').xml')

# Change the number 0.5 to 0.4 in hazard_curve-mean.xml so that it will run with the built-in parser.
try:
    lines = open(hazcurvefile, 'r').readlines()
except IOError, e:
    print
    print 'File 1 not found.'
    sys.exit(1)
lines[2] = 'xmlns="http://openquake.org/xmlns/nrml/0.4"\n'
out = open(hazcurvefile, 'w')
out.writelines(lines)
out.close()

# get annualize the curves.
curves, curlon, curlat, metadata = return_annualised_haz_curves(hazcurvefile)
imls = array(metadata['imls'])

###############################################################################
# write OQ & Frisk data to file
###############################################################################
i = 1
ii = 0
fig = plt.figure(i, figsize=(14, 10))
yhaz = 1./2475.

# set headers
jobhead = ' '.join((job.upper(), 'Mean','SA['+period+']', 'Hazard (g)'))
oqhead = '\n\n---OpenQuake---\n\nLocation,Province,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
frhead = '\n\n---2015 NBCC---\n\nLocation,Province,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
rathead = '\n\n---Ratios---\tGSCFRISK / OQ\n\nLocation,Province,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
pcdhead = '\n\n---% Difference---\tGSCFRISK / OQ\n\nLocation,Province,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
oqt = ''
frt  = ''
rat = ''
pcd = ''

# loop thru OQ curves
for lon, lat, curve in zip(curlon, curlat, curves):
    
    # loop thru FRISK curves
    for p, fh in enumerate(friskhaz):
        if fh[0] == lon and fh[1] == lat:
            # interp to frisk probs
            oqinterp = exp(interp(log(friskprobs[::-1]), log(curve[::-1]), log(imls[::-1])))[::-1]
            
            # set OQ text
            oqt += places[p] + ',' + ','.join((str('%0.4f' % x) for x in oqinterp)) + '\n'
                    
            # set frisk test
            frhaz = fh[3:-1]
            frt += places[p] + ',' + ','.join((str('%0.4f' % x) for x in frhaz)) + '\n'
            
            # get ratios
            hazrat = frhaz / oqinterp
            
            # set ratio text
            rat += places[p] + ',' + ','.join((str('%0.4f' % x) for x in hazrat)) + '\n'
            
            # calc % difference
            numer = frhaz - oqinterp
            denom = mean(vstack((frhaz, oqinterp)), axis=0)
            pcdiff = 100. * (numer / denom)
            
            # set % diff text
            pcd += places[p] + ',' + ','.join((str('%0.2f' % x) for x in pcdiff)) + '\n'
            

# combine text
outtxt = jobhead + oqhead + oqt + frhead + frt + rathead + rat + pcdhead + pcd

# write to file
csvfile =  path.join('','jobs_hazard',job, job + '_hazard_ratio_' + str(period) + '_' + job_num +'.csv')
f = open(csvfile, 'wb')
f.write(outtxt)
f.close()
            
            
            
            
            
            
            
            
            
            
            
            

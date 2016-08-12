from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
#from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
from numpy import array, exp, log, interp, loadtxt, ceil
from oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path, sep
import warnings
import sys
warnings.filterwarnings("ignore")

"""
this code is a little hokey, and only compares specialty hazard calculations run 
in both GSCFRISK and OpenQuake for single sources - file currently must be edited 
manually prior to each test
"""

###############################################################################
# settings
###############################################################################

period = '0.2'
job = 'qcss_nbcc'
job_num = '628'
fpath = 'QCSScomparison' #'CHVcomparison' #'VICMsimple' # VICMsimple, VICM_d50, VICMsimple2

'''
   397 | successful | 2016-05-11 18:38:32 | VICMbest 2015 NBCC Hazard - BEST BRANCHES
   399 | successful | 2016-05-11 19:28:32 | VICMbest 2015 NBCC Hazard - 3 GMMs
   400 | successful | 2016-05-11 19:51:35 | VICMbest 2015 NBCC Hazard - LOW GMM
   401 | successful | 2016-05-11 20:00:57 | VICMbest 2015 NBCC Hazard - HIGH GMM
   402 | successful | 2016-05-11 20:17:31 | VICMbest 2015 NBCC Hazard - 3 GMMs - REVERSE ORDER
   403 | successful | 2016-05-19 21:11:45 | VICMbest 2015 NBCC Hazard - LOW GMM
   404 | successful | 2016-05-19 22:03:13 | VICMbest 2015 NBCC Hazard - 3 GMMs - REVERSE ORDER # updated base.py and it worked with VICMsimple!   
'''

###############################################################################
# read frisk sol data
###############################################################################

#friskfile = 'OQ_'+job.upper()+'besttest'+period+'_000thperc.sol'
#friskfile = 'OQ_VICMbesttest'+period+'_000thperc.sol' # 3 GMMs
friskfile = 'OQ_CHVregion_'+period+'_20s_000thperc.sol' # full model - 20 slices
friskfile = 'OQ_QCSS_'+period+'_000thperc.sol' # full model - 20 slices
#friskfile = 'OQ_VICMbGMPE_'+period+'_000thperc.sol' # best GMM
friskpath = path.join('..','..','data','gscfrisk_comparisons',fpath,friskfile)
#friskhaz = loadtxt(friskpath, delimiter=',', skiprows=4, usecols = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
friskprobs = [0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001]

# get place names
places = []
friskhaz = []
flat = []
flon = []
lines = open(friskpath).readlines()[12:]
for line in lines:
    if len(line) > 30:
        dat = line.strip().split()
        places.append(dat[12:])
        flon.append(float(dat[0]))
        flat.append(float(dat[1]))
        friskhaz.append([float(x) for x in dat[2:12]])
        


###############################################################################
# read OQ data
###############################################################################

hazcurvefile = path.join('..','..','jobs','hazard',job,'out','hazard_curve-mean_'+job_num+'-SA('+period+').xml')
hazcurvefile = path.join('..','..','jobs','hazard',job,'out','hazard_curve-rlz-000_'+job_num+'-SA('+period+').xml')

# Change the number 0.5 to 0.4 in hazard_curve-mean.xml so that it will run with the built-in parser.
lines = open(hazcurvefile, 'r').readlines()
lines[2] = 'xmlns="http://openquake.org/xmlns/nrml/0.4"\n'
out = open(hazcurvefile, 'w')
out.writelines(lines)
out.close()

# get annualize the curves.
curves, curlon, curlat, metadata = return_annualised_haz_curves(hazcurvefile)
imls = array(metadata['imls'])

###############################################################################
# plt OQ & Frisk data
###############################################################################
i = 1
ii = 0
fig = plt.figure(i, figsize=(14, 10))
yhaz = 1./2475.

# loop thru OQ curves
for lon, lat, curve in zip(curlon, curlat, curves):
    
    # loop thru FRISK curves
    for flo, fla, fh, pl in zip(flon, flat, friskhaz, places):
        if flo == lon and fla == lat:
            ii += 1
            ax = plt.subplot(2,3,ii)
            h1 = plt.semilogy(imls, curve, 'r-', lw=2.0)
            h2 = plt.semilogy(fh, friskprobs, 'k-', lw=2.0)
            plt.title(' '.join(pl))
            plt.grid(which='both')
            plt.semilogy([0, 2.5], [yhaz, yhaz], 'k--')
            
            # get x lims from frisk haz
            thaz = exp(interp(log(1e-4), log(friskprobs[::-1]), log(fh[::-1])))

            # round to nearest 0.1
            xmax = ceil(thaz / 0.1) * 0.1
            plt.xlim([0, xmax])
            plt.ylim([1e-4, .1])
            
            if ii == 1:
                plt.legend((h1[0], h2[0]), ['OpenQuake', 'GSCFRISK'], fontsize=12)
            
            if ii == 1 or ii == 4:
                plt.ylabel('Annual Probabability of Exceedance', fontsize=14)
                
            if ii == 4 or ii == 5 or ii == 6:
                plt.xlabel(' '.join(('Mean SA('+period+' s)', 'Hazard (g)')), fontsize=14)
                
            if ii == 6:
              i += 1
              ii = 0
              plt.savefig(path.join('..','..','jobs','hazard',job,'cmp_oq_frisk_hazcurves'+str(i)+'_'+period+'_'+job_num+'.png'), format='png',bbox_inches='tight')
              fig = plt.figure(i, figsize=(14, 10))

plt.savefig(path.join('..','..','jobs','hazard',job,'cmp_oq_frisk_hazcurves'+str(i)+'_'+period+'_'+job_num+'.png'), format='png',bbox_inches='tight')
plt.show()

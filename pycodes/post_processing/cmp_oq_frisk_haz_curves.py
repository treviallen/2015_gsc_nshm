from __future__ import unicode_literals
#from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
#from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
from numpy import array, exp, log, interp, loadtxt, ceil
#from oq_tools import return_annualised_haz_curves
from hazard_tools import parse_oq_xml_poes
import matplotlib.pylab as plt
import matplotlib as mpl
mpl.style.use('classic')
from os import path, sep
import warnings, sys
import importlib
importlib.reload(sys) # for unicode chars
#sys.setdefaultencoding("latin-1")
import codecs
warnings.filterwarnings("ignore")


###############################################################################
# settings
###############################################################################
'''
period = '1.0'
job = 'secan'
job_num = '495'
'''

job = sys.argv[1]
job_num = sys.argv[2]
period = sys.argv[3]

###############################################################################
# read frisk data
###############################################################################

tmpper = period.replace('.','')
friskfile = 'NBCC2015Loc_mean_hazcurves_'+tmpper+'.csv'
friskpath = path.join('..','..','data','nbcc_mean_haz_curves',friskfile)
filecp = codecs.open(friskpath, encoding = 'latin-1')
friskhaz = loadtxt(filecp, delimiter=',', skiprows=4, usecols = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
friskprobs = [0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001]

# get place names
places = []
lines = open(friskpath, encoding = 'latin-1').readlines()[4:]
for line in lines:
    places.append(', '.join(line.strip().split(',')[-2:]))


###############################################################################
# read OQ data
###############################################################################
hazcurvefile = path.join('..','..','jobs','hazard',job,'out','hazard_curve-mean_'+job_num+'-SA('+period+').xml')

# Change the number 0.5 to 0.4 in hazard_curve-mean.xml so that it will run with the built-in parser.
#try:
"""
lines = open(hazcurvefile, 'r').readlines()
'''
except IOError, e:
    print
    print('File 1 not found.')
    sys.exit(1)
'''    
lines[2] = 'xmlns="http://openquake.org/xmlns/nrml/0.4"\n'
out = open(hazcurvefile, 'w')
out.writelines(lines)
out.close()

# get annualize the curves.
curves, curlon, curlat, metadata = return_annualised_haz_curves(hazcurvefile)
imls = array(metadata['imls'])
"""
try:
    calcDetails, calcDat = parse_oq_xml_poes(hazcurvefile)
    imls = calcDat['imls']
    curlon = calcDat['lons']
    curlat = calcDat['lats']
    curves = calcDat['annual_poes']
except:
    hazcurvefile='../../jobs/hazard/secan_collapsed_rates/secan_collapse_hazard_ratio_0.2_617_summary.csv'
    
    lines = open(hazcurvefile, encoding = 'latin-1').readlines()[1:]
    curlon = []
    curlat = []
    curves = []
    for line in lines:
        dat = line.strip().split(',')
        curlon.append(float(dat[0]))
        curlat.append(float(dat[1]))
        curves.append(array([float(x) for x in dat[4:]]))
    

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
    for p, fh in enumerate(friskhaz):
        if fh[0] == lon and fh[1] == lat:
            ii += 1
            ax = plt.subplot(2,3,ii)
            try:
                h1 = plt.semilogy(imls, curve, '-', c='darkorange', lw=2.0)
            except:
                h1 = plt.semilogy(curve, friskprobs, '-', c='darkorange', lw=2.0)
            h2 = plt.semilogy(fh[3:-1], friskprobs, '--', c='blue', lw=2.5)
            plt.title(places[p])#.encode('utf8'))
            plt.grid(which='both')
            plt.semilogy([0, 2.5], [yhaz, yhaz], 'k--')
            
            # get x lims from frisk haz
            thaz = exp(interp(log(1e-4), log(friskprobs[::-1]), log(fh[3:-1][::-1])))

            # round to nearest 0.1
            xmax = ceil(thaz / 0.1) * 0.1
            plt.xlim([0, xmax])
            plt.ylim([1e-4, .1])
            
            if ii == 1:
                plt.legend((h1[0], h2[0]), ['OpenQuake', 'GSCFRISK'], fontsize=12)
            
            if ii == 1 or ii == 4:
                plt.ylabel('Annual Probabability of Exceedance', fontsize=14)
                
            if ii == 4 or ii == 5 or ii == 6:
                plt.xlabel(' '.join(('Mean','SA['+period+']', 'Hazard (g)')), fontsize=14)
                
            if ii == 6:
              plt.savefig(path.join('..','..','jobs','hazard',job, '_'.join((job,'oq_frisk_hazcurves',period,job_num,str(i)+'.png'))), format='png',bbox_inches='tight')

              plt.savefig(path.join('..','..','jobs','hazard',job, '_'.join((job,'oq_frisk_hazcurves',period,job_num,str(i)+'.pdf'))), format='pdf',bbox_inches='tight')
              i += 1
              ii = 0
              fig = plt.figure(i, figsize=(14, 10))

if ii != 0:
    plt.savefig(path.join('..','..','jobs','hazard',job, '_'.join((job,'oq_frisk_hazcurves',period,job_num,str(i)+'.png'))), format='png',bbox_inches='tight')
    plt.savefig(path.join('..','..','jobs','hazard',job, '_'.join((job,'oq_frisk_hazcurves',period,job_num,str(i)+'.pdf'))), format='pdf',bbox_inches='tight')

plt.show()



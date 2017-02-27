from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, mean, vstack
from oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path
import warnings, sys
reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

# set NBCC probabilities
nbccprobs = array([1., 0.5, 0.2, 0.1, 0.05, 0.02,  0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001])

###############################################################################
# read config file
###############################################################################

conf_file = sys.argv[1]
leg_text  = sys.argv[2] # text to add to test curve in legend

# get paths for input files
lines = open(conf_file).readlines()
hazcurvefile1 = lines[0].split('=')[-1].strip()
hazcurvefile2 = lines[1].split('=')[-1].strip()
outputdir     = lines[2].split('=')[-1].strip()
sitelistfile  = lines[3].split('=')[-1].strip()

# get period from hazcurvefile
period = hazcurvefile1.split('(')[1].split(')')[0]

# get job numbers
job1 = hazcurvefile1.split('_')[-1].split('-')[0]
job2 = hazcurvefile2.split('_')[-1].split('-')[0]

###############################################################################
# read OQ data
###############################################################################

def get_oq_haz_curves(hazcurvefile):
    if hazcurvefile.endswith('xml'):
        # Change the number 0.5 to 0.4 in hazard_curve-mean.xml so that it will run with the built-in parser.
        lines = open(hazcurvefile, 'r').readlines()
        lines[2] = 'xmlns="http://openquake.org/xmlns/nrml/0.4"\n'
        out = open(hazcurvefile, 'w')
        out.writelines(lines)
        out.close()
    
    # get annualize the curves.
    curves, curlon, curlat, metadata = return_annualised_haz_curves(hazcurvefile)
    imls = array(metadata['imls'])
    
    return curves, curlon, curlat, metadata, imls

# get curves
curves1, curlon1, curlat1, metadata1, imls1 = get_oq_haz_curves(hazcurvefile1)
curves2, curlon2, curlat2, metadata2, imls2 = get_oq_haz_curves(hazcurvefile2)

###############################################################################
# parse site file
###############################################################################

lines = open(sitelistfile).readlines()
places = []
place_lat = []
place_lon = []

for line in lines:
    dat = line.strip().split(',')
    place_lon.append(float(dat[0]))
    place_lat.append(float(dat[1]))
    places.append(dat[2])
            
###############################################################################
# plt OQ & Frisk data
###############################################################################
i = 1
ii = 0
fig = plt.figure(i, figsize=(14, 10))
yhaz = 1./2475.

# loop thru 1st OQ curves
for lon1, lat1, curve1 in zip(curlon1, curlat1, curves1):
    
    # loop thru 2nd OQ curves
    for lon2, lat2, curve2 in zip(curlon2, curlat2, curves2):
        if lon2 == lon1 and lat2 == lat2:
            ii += 1
            ax = plt.subplot(2,3,ii)
            
            # plt haz curves
            h1 = plt.semilogy(imls1, curve1, 'r-', lw=2.0)
            h2 = plt.semilogy(imls2, curve2, '-', color='limegreen', lw=2.0)
            
            # loop thru places to get title
            for place, plon, plat in zip(places, place_lon, place_lat):
                if plat == lat1 and plon == lon1:
                    plt.title(place)#.encode('utf8'))
            
            plt.grid(which='both')
            plt.semilogy([0, 2.5], [yhaz, yhaz], 'k--')
            
            # get x lims from haz curve 2
            thaz = exp(interp(log(1e-4), log(curve2[::-1]), log(imls2[::-1])))

            # round to nearest 0.1
            xmax = ceil(thaz / 0.1) * 0.1
            plt.xlim([0, xmax])
            plt.ylim([1e-4, .1])
            plt.xlim([0, xmax])
            
            if ii == 1 or ii == 4:
                plt.ylabel('Annual Probabability of Exceedance', fontsize=14)
            """
            ###############################################################################
            # interpolate to get % diff
            ###############################################################################
            
            # interp ground motions for given probs
            
            oqinterp1 = exp(interp(log(nbccprobs[::-1]), log(curve1[::-1]), log(imls1[::-1])))[::-1]
            oqinterp2 = exp(interp(log(nbccprobs[::-1]), log(curve2[::-1]), log(imls2[::-1])))[::-1]
            
            # calculate % difference in probs & plot
            numer = oqinterp1 - oqinterp2
            denom = mean(vstack((oqinterp1, oqinterp2)), axis=0)
            pcdiff = 100. * (numer / denom)
            
            # plot twiny
            ax2 = ax.twiny()
            #ax2.set_xscale('linear')
            h3 = ax2.plot(pcdiff, nbccprobs, 'b-', lw=2.0)
            
            # plt grids manually
            ax2.plot([0, 0], [nbccprobs[0], nbccprobs[-1]], 'b--', lw=0.75)
            ax2.plot([-5, -5], [nbccprobs[0], nbccprobs[-1]], 'b--', lw=0.75)
            ax2.plot([5, 5], [nbccprobs[0], nbccprobs[-1]], 'b--', lw=0.75)
            ax2.set_xlim([-10, 10])
            
            # set col of labels
            ax2.xaxis.label.set_color('b')
            ax2.tick_params(axis='x', colors='b')
            
            if ii < 4:
                plt.xlabel('Percent Change', fontsize=13)
            """
            if ii == 1:
                plt.legend((h1[0], h2[0]), ['OQ 2015 NBCC', leg_text], fontsize=11)
            
            if ii == 4 or ii == 5 or ii == 6:
                plt.xlabel(' '.join(('Mean','SA['+period+']', 'Hazard (g)')), fontsize=13)
                                   
            if ii == 6:
              # adjust x axes
              #fig.subplots_adjust(hspace=0.2)
              
              # save
              plt.savefig(path.join(outputdir, '_'.join(('oq_hazcurves','SA('+period+')',job1,job2,str(i)+'.png'))), format='png',bbox_inches='tight')
              i += 1
              ii = 0
              fig = plt.figure(i, figsize=(14, 10))

if ii != 0:
    # adjust x axes
    #fig.subplots_adjust(hspace=0.2)
    
    plt.savefig(path.join(outputdir, '_'.join(('oq_hazcurves','SA('+period+')',job1,job2,str(i)+'.png'))), format='png',bbox_inches='tight')
plt.show()



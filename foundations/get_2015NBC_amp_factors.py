# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:50:20 2015

@author: tallen

Code generates a python dictionary of amplification factors for the NBCC2015

To access amp factors for a given PGAref, site class, type:

	> C_ampfactors[PGAref][Vs30]

	e.g. 
	
	IN > C_ampfactors['0.1g']['vs1100']
	
	OUT > array([ 0.84782554,  0.75796444,  0.69797956,  0.64357284,  0.62994067,
                0.62824066,  0.6411572 ,  0.69124828,  0.85792907,  0.67049578]) 
	
Outputs CSV file of amp factors for NBCC 2015

"""

from boore_atkinson_site_2008 import boore_atkinson_siteamp
from atkinson_boore_site_2006 import atkinson_boore_siteamp
from numpy import array, arange, where, log, exp, interp, hstack

# periods for NBCC2015 (PGA = 0; PGV = -1)     
T = [0.05, 0.1, 0.2, 0.3, 0.5, 1., 2., 3., 5., 7.5, 10., 0., -1.] # in s
#T = [0.1, 0.2, 0.3, 0.5, 1., 2., 5., 10., 0., -1.] # in s - NBCC


# Set Vs30 array
vs30 = [115, 250, 450, 1100, 1600] # in m/s 

# frist get B/C pga4nl
target_pga4nl_C = array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5]) # in g

# get B/C linear ground-motions that provide target_pga4nl_C values - AW rounding here led to minor differences in final factors
refVs30 = 450.
refT    = 0.0 # i.e. PGA
refPGA  = 0.1 # in g
pga4nl_BC = target_pga4nl_C / boore_atkinson_siteamp(refVs30, refT, refPGA)[0] # B/C to C conversion factor for PGA of 1.208 based on BA08

#pga4nl_BC = target_pga4nl_C / 1.  # for testing purposes - needed to match Gail's BA08 C->B/C factors

##########################################################################
'''
Get BA08 amp factors for site class BC - this loop is not required for 
NBCC2015 site class tables, but kept for testing purposes
'''

# set empty dictionary
BC_ampfactors = {}

# loop thru pga4nl
for i, pga in enumerate(pga4nl_BC):
    # set dictionary and loop thru velocity
    vdict = {}
    # loop thru Vs30
    for v in vs30:
        # loop thru periods
        amp = []
        for t in T:
            # get amp factor for each period and append to amp
            amp.append(boore_atkinson_siteamp(v, t, pga)[0])
        
        # add amp factors to velocity dictionary
        vdict['vs'+str(v)] = array(amp)
    
    # add velocity dictionary to pga4nl dictionary
    BC_ampfactors[str(target_pga4nl_C[i])+'g'] = vdict
    
BC_ampfactors['periods'] = array(T)

##########################################################################
'''
Re-cast amp factors in terms of site class C 
This gets appropriate amp factors for classes C-E
'''
# set empty dictionary
C_ampfactors  = {}

# loop thru pga4nl
for i, pga in enumerate(pga4nl_BC):
    # set dictionary and loop thru velocity
    vdict = {}
    for v in vs30:
        # loop thru periods
        amp = []
        for t in T:
            # get amp factor for each period relative to C and append to amp
            amp.append(boore_atkinson_siteamp(v, t, pga)[0] / 
                       boore_atkinson_siteamp(450., t, pga)[0])
        
        # add amp factors to velocity dictionary
        vdict['vs'+str(v)] = array(amp)
    
    # add velocity dictionary to amp factor dictionary
    C_ampfactors[str(target_pga4nl_C[i])+'g'] = vdict

# add periods to amp factor dictionary 
C_ampfactors['periods'] = array(T)

##########################################################################
'''
For sites faster than 760 m/s, use GA's proposed interpolation to obtain 
factors for A-B site classes (1600 & 1100 m/s)
'''

# BA08 factors to convert C (vs30 = 450) to BC (vs30 = 760)
GA_T = [0.1, 0.2, 0.3, 0.5, 1., 2., 5., 10., 0., -1.] # in s
BA08_CtoBC = [0.88, 0.85, 0.79, 0.73, 0.69, 0.68, 0.68, 0.71, 0.83, 0.73] # from Gail's spreadsheet - not used!
# interp from GA periods to NBCC periods
BA08_CtoBC = hstack((exp(interp(log(T[0:-2]), log(GA_T[0:-2]), log(BA08_CtoBC[0:-2]))), \
                    BA08_CtoBC[-2:]))

'''
# BA08 factors to convert C (vs30 = 450) to BC (vs30 = 760) using pga4nl value for BC that gives 0.1g on C
BA08_CtoBC = []
for t in T:
    BA08_CtoBC.append(boore_atkinson_siteamp(760., t, pga4nl_BC[0]) \
                      / boore_atkinson_siteamp(450., t, pga4nl_BC[0])) # amp factors determined using PGAref[C] = 0.1 g
BA08_CtoBC = array(BA08_CtoBC)
'''

# AB06 factors to convert C to vs30 = 2000 m/s
AB06_Cto2000 = [0.82, 0.65, 0.58, 0.53, 0.54, 0.55, 0.59, 0.66, 0.93, 0.59]
# interp from GA periods to NBCC periods
AB06_Cto2000 = hstack((exp(interp(log(T[0:-2]), log(GA_T[0:-2]), log(AB06_Cto2000[0:-2]))), \
                    AB06_Cto2000[-2:]))

# reference Vs30
ref_vs30 = [760., 2000.]

# target Vs30
tar_vs30 = [1100, 1600]

# log-log interpolate between amp factors - loop thru periods
for i in range(0, len(T)):
    # set amplitudes to interpolate between
    interp_amps = [BA08_CtoBC[i], AB06_Cto2000[i]]
    
    # do log-log interpolation
    interp_factors = exp(interp(log(tar_vs30), log(ref_vs30), log(interp_amps)))
                          
    # supplant newly interpolated factors for classes A-B - factors the same for all PGAref
    for pga in target_pga4nl_C:
        C_ampfactors[str(pga)+'g']['vs'+str(tar_vs30[0])][i] = interp_factors[0]
        C_ampfactors[str(pga)+'g']['vs'+str(tar_vs30[1])][i] = interp_factors[1]

##########################################################################
'''
Get CSV text for output to match NBCC tables
'''

# set params
ascii_no = range(65, 70) # ascii chars A-E
site_class = [chr(x) for x in ascii_no] # convert ascii number to char
rev_vs30 = vs30[::-1] # reverse vs30 order for tables
PGAref_header = ','+','.join(['PGAref='+str(x)+'g' for x in target_pga4nl_C])+'\n'
NBCCtxt = ''

# loop thru period
for i, t in enumerate(T):
    # get period header
    T_header = ','.join(('Site Class', 'Values of F(T) for T = '+str(t)+' s\n'))
    
    # add headers to NBCC text
    NBCCtxt += T_header + PGAref_header
    
    # add amp factors for each site class
    for j, sc in enumerate(site_class):
        # loop thru PGAref and add amp factor
        site_txt = ''
        for pga in target_pga4nl_C:
            # round to 2 decimal points for output
            site_txt += ',' + str('%0.2f' % C_ampfactors[str(pga)+'g']['vs'+str(rev_vs30[j])][i])
        
        # add site class txt to NBCC text
        NBCCtxt += sc + site_txt + '\n'
    NBCCtxt += '\n'

# write NBCC factors to file
f = open('NBCC2015_ampfactors.csv', 'wb')
f.write(NBCCtxt)
f.close()

##########################################################################
'''
plt amp factors
'''

import matplotlib.pyplot as plt
import matplotlib as mpl

pltpga = [0.1, 0.4]

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14 

numcols = 5
cmap = plt.cm.get_cmap('rainbow', numcols)
cs = (cmap(arange(numcols)))

# remove yellow
#idx=array([0, 1, 3, 4, 5])
#cs = cs[idx]

legtxt = ['Site Class A', 'Site Class B', 'Site Class C', 'Site Class D', 'Site Class E', 'F(T)', 'Fa & Fv']
figure = plt.figure(1,figsize=(19,8))

# Boarcherdt 1994 for 0.1 and 0.4 g
Fa01 = [0.6, 0.7, 1.0, 1.5, 1.5]
Fa04 = [1.1, 1.0, 1.0, 0.9, 0.9]
Fv01 = [0.4, 0.6, 1.0, 2.0, 2.0]
Fv04 = [0.6, 0.7, 1.0, 1.6, 1.6]
FaT  = [0.1, 0.5]
FvT  = [0.4, 2.0]


# Finn & Whiteman for 0.1 and 0.4 g - Relative to site Sa(0.2 s)
Fa01 = array([0.7, 0.8, 1.0, 1.3, 2.1])  
Fa04 = array([1.1, 1.0, 1.0, 0.9, 0.9])  
Fv01 = array([0.4, 0.6, 1.0, 2.0, 2.0])  
Fv04 = array([0.6, 0.7, 1.0, 1.6, 1.6])  
FaT = [0.1, 0.5]
FvT = [0.4, 2.0] 


for i, pga in enumerate(pltpga):
    ax = plt.subplot(1, 2, i+1)
    plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    # loop thru vs
    for j, v in enumerate(rev_vs30):
        pltamps = C_ampfactors[str(pga)+'g']['vs'+str(v)][:-2]
        plt.semilogx(T[:-2], pltamps, '-', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])        
        	
    plt.ylim([0.5, 3.])
    plt.xlim([0.05, 10.])
    plt.title('PGAref = ' + str(pga) + ' g', fontsize=18)
    plt.xlabel('Period (s)', fontsize=16)
    plt.ylabel('Amplification Factor (Relative to C)', fontsize=16)
    
    '''
    if i == 0:
        plt.legend(legtxt ,loc=2, fontsize=14)
    else:
        h1 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k-', lw=3.)
        h2 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k--', lw=3.)
        plt.legend((h1[0], h2[0]),('F(T)', 'Fa & Fv') ,loc=2, fontsize=14)
    '''
    # at John's request
    plt.semilogx([0.1, 1.],[9999., 9999.], 'k-', lw=4.)
    plt.semilogx([0.1, 1.],[9999., 9999.], 'k--', lw=4.)
    plt.legend(legtxt ,loc=2, fontsize=13)
        
    for j, v in enumerate(rev_vs30):	
        # plt Fa, Fv
        if i == 0:
            plt.semilogx(FaT, [Fa01[j], Fa01[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            plt.semilogx(FvT, [Fv01[j], Fv01[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
        else:                                              
            plt.semilogx(FaT, [Fa04[j], Fa04[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            plt.semilogx(FvT, [Fv04[j], Fv04[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            
    # set new xticls
    ax.set_xticks([0.1, 0.3, 1., 3., 10])
    xlabels = ['0.1', '0.3', '1.0', '3.0', '10']
    ax.set_xticklabels(xlabels)
    plt.ylim([0.5, 4.1])
        
plt.savefig('2015NBCC_F(T)_factors.png', format='png', bbox_inches='tight', dpi=150)

plt.show()
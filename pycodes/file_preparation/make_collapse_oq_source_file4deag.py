#####################################################################################
# reads pkl file of 2015NBCC hazard model and outputs as OpenQuake source files 
# for sources within a given distance of the site, as well as the source logic tree
#
# NOTE: lines 86-89 should be modified for each site of interest

# Usage: python make_collapse_oq_source_file_deag.py <path to pkl file> <include src weights>
#            e.g. python make_collapse_oq_source_file4deag.py ../../data/pkl_files/SECan_H2_src_model.pkl True
#####################################################################################

def make_collapse_occurrence_text(m, binwid):
    from numpy import zeros    
    
    bval_wt    = [0.68, 0.16, 0.16]
    max_mag_wt = [0.60, 0.30, 0.10]
    wtd_list  = []
    maglen = 0

    for i, bwt in enumerate(bval_wt):
        # get effective rate
        effN0 = m['src_N0'][i] * m['src_weight']
        
        for j, mwt in enumerate(max_mag_wt):
            
            betacurve, mrange = get_oq_incrementalMFD(m['src_beta'][i], effN0, \
                                                      m['min_mag'], m['max_mag'][j], \
                                                      binwid)
            wtd_list.append(betacurve * bwt * mwt)  
            
            # get max length of arrays
            if len(betacurve) > maglen:
                maglen = len(betacurve)

    # sum rates for each branch
    wtd_rates = zeros(maglen)
    for rates in wtd_list:
        # go element by element
        for r, rate in enumerate(rates):
            wtd_rates[r] += rate
                
    # convert cummulative rates to annual occurrence rates
    occ_rates = []
    for b in range(0, len(wtd_rates[0:-1])):
        occ_rates.append(wtd_rates[b] - wtd_rates[b+1])
    occ_rates.append(wtd_rates[-1])
    
    # make text object                        
    octxt = str('%0.5e' % occ_rates[0])
    for bc in occ_rates[1:]:
        octxt += ' ' + str('%0.5e' % bc)
        
    return octxt


'''
start main code here
'''

#import datetime as dt
from sys import argv
import pickle
from oq_tools import beta2bval, get_oq_incrementalMFD, get_line_parallels, distance
from numpy import array, log10, max, min, tan, radians, unique, isinf, isnan, concatenate
from os import path

########################################################################################

# this section needs to be edited
sitename = 'QCFsite'
sitelat = 53.45
sitelon = -132.80
maxgmpedist = 200.

########################################################################################

inputpkl  = argv[1]
multimods = argv[2] # for setting weights of alternative models (True or False)
# site specific - ignore #outbase   = argv[3] # folder for github sources to be included in source_model_logic_tree.xml

# read pickle
pklfile = open(inputpkl, 'rb')
model = pickle.load(pklfile)

# set big bbox params
bbmaxlon = -180
bbmaxlat = -90
bbminlon = 180
bbminlat = 90

# Write 9 model files
betalist = ['bb','bl','bu']
maglist = ['mb','ml','mu']
srcxmls = []

# make xml header
header = '<?xml version="1.0" encoding="utf-8"?>\n'
header += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
header += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'

# set wieghts
bval_wt    = [0.68, 0.16, 0.16]
max_mag_wt = [0.60, 0.30, 0.10]
branch_wt = []

# loop through seiscomp xml template and make sourcefile
outbase = model[0]['src_reg'].replace(' ','_')
newxml = header + '    <sourceModel name="'+outbase+'">\n\n'

# get src codes and rename if duplicated
codes = []
for m in model:
    codes.append(m['src_code'])
ucodes = unique(codes)

for m in model:
    
    # get tectonic region
    if m['gmpe'] == 'Wcrust_med_clC.txt':
        tectonic_reg = 'Active Shallow Crust'
    elif m['gmpe'] == 'WinterfaceCombo_medclC.txt':
        tectonic_reg = 'Subduction Interface'
    elif m['gmpe'] == 'ENA_med_clC.txt':
        tectonic_reg = 'Stable Shallow Crust'
    elif m['gmpe'] == 'WcrustFRjb_med_clC.txt':
        tectonic_reg = 'Active Shallow Fault'
    elif m['gmpe'] == 'Woffshore_med_clC.txt':
        tectonic_reg = 'Active Shallow Offshore'
    elif m['gmpe'] == 'WinslabD30_med_clC.txt':
        tectonic_reg = 'Subduction IntraSlab30'
    elif m['gmpe'] == 'WinslabD50_med_clC.txt':
        tectonic_reg = 'Subduction IntraSlab50'
        
    # check to see if source is in valid distance range
    writesrc = False
    srclola = m['src_shape']
    for lo, la in srclola:
        rngkm, az, baz = distance(sitelat, sitelon, la, lo)
        
        # check if in distance rng
        if rngkm <= maxgmpedist:
            writesrc = True
    
    #######################################################################
    # write area sources
    #######################################################################
    if m['src_type'] == 'area' and writesrc == True:
        print m['src_type']
        
        # rename source code if "." exists
        m['src_code'].replace('.', '')
        
        # ignore SCC sources
        #if m['src_code'].startswith('SCC') == False:
        #if tectonic_reg != 'Stable Shallow Crust':
    
        newxml += '        <areaSource id="'+m['src_code']+'" name="'+\
                   m['src_name']+'" tectonicRegion="'+tectonic_reg+'">\n'
        
        newxml += '            <areaGeometry>\n'
        newxml += '                <gml:Polygon>\n'
        newxml += '                    <gml:exterior>\n'
        newxml += '                        <gml:LinearRing>\n'
        newxml += '                            <gml:posList>\n'
                        
        # get polygon text
        polytxt = ''
        for xy in m['src_shape'][:-1]: # no need to close poly
            polytxt = polytxt + '                                ' + str("%0.4f" % xy[0]) \
                              + ' ' + str("%0.4f" % xy[1]) + '\n'
        newxml += polytxt
        
        newxml += '                            </gml:posList>\n'
        newxml += '                        </gml:LinearRing>\n'
        newxml += '                    </gml:exterior>\n'
        newxml += '                </gml:Polygon>\n'
        
        ###################################################################
        # print model bbox of model
        
        # this is not required for the nrml files, but useful for setting up job.ini files
        
        buff = 0.1
        maxlon = max(m['src_shape'][:,0])+buff
        minlon = min(m['src_shape'][:,0])-buff
        maxlat = max(m['src_shape'][:,1])+buff
        minlat = min(m['src_shape'][:,1])-buff

        # get big bbox
        if maxlon > bbmaxlon: bbmaxlon = maxlon
        if minlon < bbminlon: bbminlon = minlon
        if maxlat > bbmaxlat: bbmaxlat = maxlat
        if minlat < bbminlat: bbminlat = minlat

        print m['src_code'], minlon, minlat, ',', minlon, maxlat, ',', maxlon, maxlat, ',', maxlon, minlat
        ###################################################################

        # set depth distribution
        if min(m['src_dep']) != max(m['src_dep']):
            newxml += '                <upperSeismoDepth>'+str("%0.1f" % min(m['src_dep']))+'</upperSeismoDepth>\n'
            newxml += '                <lowerSeismoDepth>'+str("%0.1f" % max(m['src_dep']))+'</lowerSeismoDepth>\n'
        else:
            newxml += '                <upperSeismoDepth>'+str("%0.1f" % (min(m['src_dep'])-10))+'</upperSeismoDepth>\n'
            newxml += '                <lowerSeismoDepth>'+str("%0.1f" % (min(m['src_dep'])+10))+'</lowerSeismoDepth>\n'
            
        newxml += '            </areaGeometry>\n'
        newxml += '            <magScaleRel>PointMSR</magScaleRel>\n'
        newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
        
        # get weighted rates
        binwid = 0.1
        octxt = make_collapse_occurrence_text(m, binwid)
        
        newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
        newxml += '                <occurRates>'+octxt+'</occurRates>\n'
        newxml += '            </incrementalMFD>\n'
        
        # set nodal planes
        newxml += '            <nodalPlaneDist>\n'
        newxml += '                <nodalPlane probability="1.0" strike="0.0" dip="90.0" rake="0.0" />\n' # not necessary for replicating 2015NBCC
        """
        newxml += '                <nodalPlane probability="0.125" strike="0.0" dip="90.0" rake="0.0" />\n'
        newxml += '                <nodalPlane probability="0.125" strike="45.0" dip="90.0" rake="0.0" />\n'
        newxml += '                <nodalPlane probability="0.125" strike="90.0" dip="90.0" rake="0.0" />\n'
        newxml += '                <nodalPlane probability="0.125" strike="135.0" dip="90.0" rake="0.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="0.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="45.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="90.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="135.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="180.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="225.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="270.0" dip="30.0" rake="90.0" />\n'
        newxml += '                <nodalPlane probability="0.0625" strike="315.0" dip="30.0" rake="90.0" />\n'
        """
        newxml += '            </nodalPlaneDist>\n'
        

        # set hypo depth
        newxml += '            <hypoDepthDist>\n'
        newxml += '                <hypoDepth probability="0.50" depth="'+str("%0.1f" % m['src_dep'][0])+'"/>\n' \
                 +'                <hypoDepth probability="0.25" depth="'+str("%0.1f" % m['src_dep'][1])+'"/>\n' \
                 +'                <hypoDepth probability="0.25" depth="'+str("%0.1f" % m['src_dep'][2])+'"/>\n'
        newxml += '            </hypoDepthDist>\n'
        newxml += '        </areaSource>\n\n'

    #######################################################################
    # now make fault sources
    #######################################################################
    elif m['src_type'] == 'fault' and writesrc == True:
        
        # rename source code if "." exists
        m['src_code'].replace('.', '')
        
        #print m['src_code'], log10(m['src_N0'][i])
        if isinf(log10(m['src_N0'][0])) == False:
            ###################################################################
            # do complex faults
            ###################################################################
            if m['fault_dip'][0] != m['fault_dip'][1]: # old
            #if m['fault_dip'][0] >= 0: # catches all faults
                #if m['fault_dip'][0] > 0:
                if m['src_code'].startswith('CASCADIA'):
                    src_code = 'CIS'
                else:
                    src_code = m['src_code']
                    
                # id subcript
                idsub = str("%0.1f" % beta2bval(m['src_beta'][0]))
                idsub = idsub.replace(".", "")
                
                newxml += '        <complexFaultSource id="'+src_code+idsub+'" name="'+\
                           m['src_name']+'" tectonicRegion="'+tectonic_reg+'">\n'
                newxml += '            <complexFaultGeometry>\n'
                newxml += '                <faultTopEdge>\n'
                newxml += '                    <gml:LineString>\n'
                newxml += '                        <gml:posList>\n'
            
                # calculate lat lons from surface projection
                # get upper h-dist
                upperhdist = m['src_dep'][0] / tan(radians(m['fault_dip'][0]))
                upperxy = get_line_parallels(m['src_shape'], upperhdist)[0]
            
                # make upper text
                xytxt = ''
                for xy in upperxy:
                    xytxt += '                            ' + \
                             ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1]), str(m['src_dep'][0])))+'\n'
                newxml += xytxt
                newxml += '                        </gml:posList>\n'
                newxml += '                    </gml:LineString>\n'
                newxml += '                </faultTopEdge>\n'
                newxml += '                <intermediateEdge>\n'
                newxml += '                    <gml:LineString>\n'
                newxml += '                        <gml:posList>\n'
            
                
                # calculate lat lons from upper edge
                # get intermediate h-dist
                interhdist = (m['src_dep'][1] - m['src_dep'][0]) / tan(radians(m['fault_dip'][0]))
                interxy = get_line_parallels(upperxy, interhdist)[0]
            
                # make intermediate text
                xytxt = ''
                for xy in interxy:
                    xytxt += '                            ' + \
                             ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1]), str(m['src_dep'][1])))+'\n'
                newxml += xytxt
                newxml += '                        </gml:posList>\n'
                newxml += '                    </gml:LineString>\n'
                newxml += '                </intermediateEdge>\n'
                newxml += '                <faultBottomEdge>\n'
                newxml += '                    <gml:LineString>\n'
                newxml += '                        <gml:posList>\n'
            
                # calculate lat lons from intermediate edge
                # get bottom h-dist
                bottomhdist = (m['src_dep'][2] - m['src_dep'][1]) / tan(radians(m['fault_dip'][1]))
                bottomxy = get_line_parallels(interxy, bottomhdist)[0]
            
                # make bottom text
                xytxt = ''
                for xy in bottomxy:
                    xytxt += '                            ' + \
                             ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1]), str(m['src_dep'][2])))+'\n'
                newxml += xytxt
                newxml += '                        </gml:posList>\n'
                newxml += '                    </gml:LineString>\n'
                newxml += '                </faultBottomEdge>\n'
                newxml += '            </complexFaultGeometry>\n'
                
                '''
                # get fault area scaling model
                '''
                #src_code = m['src_code']
                if src_code.startswith('CIS'):
                    newxml += '            <magScaleRel>GSCCascadia</magScaleRel>\n'
                elif src_code.startswith('WIN'):
                    newxml += '            <magScaleRel>GSCOffshoreThrustsWIN</magScaleRel>\n'
                elif src_code.startswith('HGT'):
                    newxml += '            <magScaleRel>GSCOffshoreThrustsHGT</magScaleRel>\n'
                elif src_code.startswith('QCSS') or src_code.startswith('FWF'):
                    newxml += '            <magScaleRel>WC1994_QCSS</magScaleRel>\n'
                elif src_code.startswith('EISO'):
                    newxml += '            <magScaleRel>GSCEISO</magScaleRel>\n'
                elif src_code.startswith('EISB'):
                    newxml += '            <magScaleRel>GSCEISB</magScaleRel>\n'
                elif src_code.startswith('EISI'):
                    newxml += '            <magScaleRel>GSCEISI</magScaleRel>\n'
                else:
                    newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
                
                newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
            
                '''
                # now get appropriate MFD
                '''
                # do incremental MFD
                if m['src_beta'][0] > -99:
                    # set incremental recurrence pars
                    octxt = make_collapse_occurrence_text(m, binwid)
            
                    # make text
                    newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                    newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                    newxml += '            </incrementalMFD>\n'
            
                if m['fault_dip'][0] != 90.:
                    newxml += '            <rake>90.0</rake>\n'
                else:
                    newxml += '            <rake>0.0</rake>\n'
            
                newxml += '        </complexFaultSource>\n\n'

            ###################################################################
            # else do simple fault
            ###################################################################
            elif m['fault_dip'][0] == m['fault_dip'][1]:
                
                # id subcript
                idsub = str("%0.1f" % beta2bval(m['src_beta'][0]))
                idsub = idsub.replace(".", "")
                
                newxml += '        <simpleFaultSource id="'+m['src_code']+idsub+'" name="'+\
                                     m['src_name']+'" tectonicRegion="'+tectonic_reg+'">\n'
                newxml += '            <simpleFaultGeometry>\n'
                newxml += '                <gml:LineString>\n'
                newxml += '                    <gml:posList>\n'
            
                # simple fauls use surface projection!
                '''
                # calculate lat lons from surface projection
                # get upper h-dist
                upperhdist = m['src_dep'][0] / tan(radians(m['fault_dip'][0]))
                upperxy = get_line_parallels(m['src_shape'], upperhdist)[0]
                '''
                
                xytxt = ''
                for xy in m['src_shape']:
                    xytxt += '                            ' + \
                             ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1])))+'\n'
                newxml += xytxt
                
                newxml += '                    </gml:posList>\n'
                newxml += '                </gml:LineString>\n'
                newxml += '                <dip>'+str(m['fault_dip'][0])+'</dip>\n'
                newxml += '                <upperSeismoDepth>'+str(m['src_dep'][0])+'</upperSeismoDepth>\n'
                newxml += '                <lowerSeismoDepth>'+str(m['src_dep'][-1])+'</lowerSeismoDepth>\n'
                newxml += '            </simpleFaultGeometry>\n'
                
                '''
                # get fault area scaling model
                '''
                src_code = m['src_code']
                if src_code == 'CIS':
                    newxml += '            <magScaleRel>GSCCascadia</magScaleRel>\n'
                elif src_code.startswith('WIN'):
                    newxml += '            <magScaleRel>GSCOffshoreThrustsWIN</magScaleRel>\n'
                elif src_code.startswith('HGT'):
                    newxml += '            <magScaleRel>GSCOffshoreThrustsHGT</magScaleRel>\n'
                elif src_code.startswith('QCSS') or src_code.startswith('FWF'):
                    newxml += '            <magScaleRel>WC1994_QCSS</magScaleRel>\n'
                elif src_code.startswith('EISO'):
                    newxml += '            <magScaleRel>GSCEISO</magScaleRel>\n'
                elif src_code.startswith('EISB'):
                    newxml += '            <magScaleRel>GSCEISB</magScaleRel>\n'
                elif src_code.startswith('EISI'):
                    newxml += '            <magScaleRel>GSCEISI</magScaleRel>\n'
                else:
                    newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
                
                newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
                #newxml += '            <ruptAspectRatio>2.0</ruptAspectRatio>\n'
                '''
                # now get appropriate MFD
                '''
                # do incremental MFD
                if m['src_beta'][0] > -99:
                    # set incremental recurrence pars
                    octxt = make_collapse_occurrence_text(m, binwid)
            
                    # make text
                    newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                    newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                    newxml += '            </incrementalMFD>\n'
            
                if m['fault_dip'][0] != 90.:
                    newxml += '            <rake>90.0</rake>\n'
                else:
                    newxml += '            <rake>0.0</rake>\n'
            
                newxml += '        </simpleFaultSource>\n\n'
            
# finish nrml
newxml += '    </sourceModel>\n'
newxml += '</nrml>'

# write Big BBOX
print '\nBBOX:', bbminlon, bbminlat, ',', bbminlon, bbmaxlat, ',', bbmaxlon, bbmaxlat, ',', bbmaxlon, bbminlat

# write new data to file
outxml = path.join('temp_source_files', '_'.join((sitename,outbase,'deag.xml')))
f = open(outxml,'w')
f.write(newxml)
f.close() 

srcxmls.append('_'.join((sitename,outbase,'deag.xml')))

######################################################################
# now that the source file have been written, make the logic tree file
# if multimodel - adjust weights
if multimods == 'True':
    branch_wt = array(branch_wt)
    branch_wt *= m['src_reg_wt']
    print 'Branch Weights: ', m['src_reg_wt']

newxml = '<?xml version="1.0" encoding="UTF-8"?>\n'
newxml += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
newxml += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'
newxml += '    <logicTree logicTreeID="lt1">\n'
newxml += '        <logicTreeBranchingLevel branchingLevelID="bl1">\n'
newxml += '            <logicTreeBranchSet uncertaintyType="sourceModel"\n' \
          '                                branchSetID="bs1">\n\n'

# make branches
for i, branch in enumerate(srcxmls):
    newxml += '                <logicTreeBranch branchID="b' + str(i+1) + '">\n'
    newxml += '                    <uncertaintyModel>'+branch+'</uncertaintyModel>\n'
    newxml += '                    <uncertaintyWeight>1.0</uncertaintyWeight>\n'
    newxml += '                </logicTreeBranch>\n\n'

newxml += '            </logicTreeBranchSet>\n'
newxml += '        </logicTreeBranchingLevel>\n'
newxml += '    </logicTree>\n'
newxml += '</nrml>'
        
# write logic tree to file
outxml = path.join('temp_source_files', ''.join((sitename, '_deag_source_model_logic_tree.xml')))
f = open(outxml,'w')
f.write(newxml)
f.close()

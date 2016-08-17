#####################################################################################
# reads pkl file of 2015NBCC hazard model and outputs as OpenQuake source files 
# for best, upper & lower models as well as the source logic tree

# Usage: python make_openquake_source_file.py <path to pkl file> <output base name>
#            e.g. python make_openquake_source_file.py SWCan_T3EclC1.pkl swcan
#####################################################################################

def make_incremental_occurrence_text(beta, effN0, mmin, mmax, binwid):
    betacurve, mrange = get_oq_incrementalMFD(beta, effN0, mmin, mmax, binwid)
    
    # convert cummulative rates to annual occurrence rates
    occ_rates = []
    for b in range(0, len(betacurve[0:-1])):
        occ_rates.append(betacurve[b] - betacurve[b+1])
    occ_rates.append(betacurve[-1])
    
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
#import shapefile
import pickle
from oq_tools import beta2bval, get_oq_incrementalMFD, get_line_parallels
from numpy import array, log10, max, min, tan, radians, unique, isinf, isnan, concatenate
from os import path

inputpkl = argv[1]
multimods = argv[2] # for setting weights of alternative models (True or False)
src_folder = argv[3] # folder for github sources to be included in source_model_logic_tree.xml

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

for i, bl in enumerate(betalist):
    for j, ml in enumerate(maglist):
        # loop through seiscomp xml template and make sourcefile
        outbase = model[0]['src_reg'].replace(' ','_')

        # fix ugly naming convention
        if outbase == '20130408_EA_H_Model_v5':
            outbase = 'E_ARCTIC_H'
        elif outbase == '20130408_EA_R_Model_v3':
            outbase = 'E_ARCTIC_R'
        elif outbase == '20131002_WA_H_Model_v7':
            outbase = 'W_ARCTIC'
        
        newxml = header + '    <sourceModel name="'+outbase+'_'+betalist[i]+'_'+maglist[j]+'">\n\n'
        
        # get src codes and rename if duplicated
        codes = []
        for m in model:
            codes.append(m['src_code'])
        ucodes = unique(codes)
        
        # get wieght for branch
        branch_wt.append(bval_wt[i] * max_mag_wt[j])
    
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
            
            #######################################################################
            # write area sources
            #######################################################################
            if m['src_type'] == 'area':
                
                # rename source code if "." exists
                m['src_code'].replace('.', '')
                
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
                
                # set incremental recurrence pars
                binwid = 0.1
                effN0 = m['src_N0'][i] * m['src_weight']
                
                octxt = make_incremental_occurrence_text(m['src_beta'][i], effN0, \
                                                         m['min_mag'], m['max_mag'][j], \
                                                         binwid)
                                                         
                newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                newxml += '            </incrementalMFD>\n'
                
                """
                # set GR recurrence pars
                tmpN0   = m['src_N0'][i]
                tmpbeta = m['src_beta'][i]
                tmpmmax = m['max_mag'][j]
                grtxt = ''.join(('            <truncGutenbergRichterMFD aValue="', \
                                str("%0.4f" % log10(tmpN0)),'" bValue="', \
                                str("%0.4f" % beta2bval(tmpbeta)),'" minMag="', \
                                str("%0.2f" % m['min_mag']),'" maxMag="', \
                                str("%0.2f" % tmpmmax),'"/>\n'))
                                    
                newxml += grtxt
                """
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
            elif m['src_type'] == 'fault':
                
                # rename source code if "." exists
                m['src_code'].replace('.', '')
                
                print m['src_code'], log10(m['src_N0'][i])
                if isinf(log10(m['src_N0'][i])) == False:
                    ###################################################################
                    # do complex faults
                    ###################################################################
                    if m['fault_dip'][0] != m['fault_dip'][1]:
                    #if m['fault_dip'][0] >= 0: # catches all faults
                        #if m['fault_dip'][0] > 0:
                        if m['src_code'].startswith('CASCADIA'):
                            src_code = 'CIS'
                        else:
                            src_code = m['src_code']
                            
                        # id subcript
                        idsub = str("%0.1f" % beta2bval(m['src_beta'][i]))
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
                        if m['src_beta'][i] > -99:
                            # adjust N0 value to account for weighting of fault sources
                            effN0 = m['src_N0'][i] * m['src_weight']
                            
                            octxt = make_incremental_occurrence_text(m['src_beta'][i], effN0, \
                                                                     m['min_mag'], m['max_mag'][j], \
                                                                     binwid)
                    
                            # make text
                            newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                            newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                            newxml += '            </incrementalMFD>\n'
                    
                        # do GR MFD
                        elif m['src_beta'][i] >= 0:
                            # set recurrence pars
                            effN0   = m['src_N0'][i] * m['src_weight']
                            tmpbeta = m['src_beta'][i]
                            tmpmmax = m['max_mag'][j]
                            if tmpbeta > 9999.:
                                mfdtype = 'YoungsCoppersmith1985MFD'
                            else:
                                mfdtype = 'truncGutenbergRichterMFD'
                                newxml += ''.join(('            <',mfdtype,' aValue="', \
                                                str("%0.4f" % log10(effN0)),'" bValue="', \
                                                str("%0.4f" % beta2bval(tmpbeta)),'" minMag="', \
                                                str("%0.2f" % m['min_mag']),'" maxMag="', \
                                                str("%0.2f" % tmpmmax),'"/>\n'))
                    
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
                        idsub = str("%0.1f" % beta2bval(m['src_beta'][i]))
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
                        if m['src_beta'][i] > -99:
                            # adjust N0 value to account for weighting of fault sources
                            effN0 = m['src_N0'][i] * m['src_weight']
                            
                            octxt = make_incremental_occurrence_text(m['src_beta'][i], effN0, \
                                                                     m['min_mag'], m['max_mag'][j], \
                                                                     binwid)
                    
                            # make text
                            newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                            newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                            newxml += '            </incrementalMFD>\n'
                    
                        # do GR MFD
                        elif m['src_beta'][i] >= 0:
                            # set recurrence pars
                            effN0   = m['src_N0'][i] * m['src_weight']
                            tmpbeta = m['src_beta'][i]
                            tmpmmax = m['max_mag'][j]
                            mfdtype = 'truncGutenbergRichterMFD'
                            newxml += ''.join(('            <',mfdtype,' aValue="', \
                                               str("%0.4f" % log10(effN0)),'" bValue="', \
                                               str("%0.4f" % beta2bval(tmpbeta)),'" minMag="', \
                                               str("%0.2f" % m['min_mag']),'" maxMag="', \
                                               str("%0.2f" % tmpmmax),'"/>\n'))
                                                
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
        outxml = path.join('temp_source_files', '_'.join((outbase,bl,ml+'.xml')))

        f = open(outxml,'w')
        f.write(newxml)
        f.close() 
        
        srcxmls.append(outxml)
        
######################################################################
# now that the source file have been written, make the logic tree file
#if multimodel - adjust weights:
if multimods == True:
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
    logictreepath = src_folder + '/' + path.split(branch)[-1]
    newxml += '                <logicTreeBranch branchID="b' + str(i+1) + '">\n'
    newxml += '                    <uncertaintyModel>'+logictreepath+'</uncertaintyModel>\n'
    newxml += '                    <uncertaintyWeight>'+str('%0.3f' % branch_wt[i])+'</uncertaintyWeight>\n'
    newxml += '                </logicTreeBranch>\n\n'

newxml += '            </logicTreeBranchSet>\n'
newxml += '        </logicTreeBranchingLevel>\n'
newxml += '    </logicTree>\n'
newxml += '</nrml>'
        
# write logic tree to file
outxml = path.join('temp_source_files', ''.join((outbase, '_source_model_logic_tree.xml')))
f = open(outxml,'w')
f.write(newxml)
f.close()

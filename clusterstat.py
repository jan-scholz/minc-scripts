#!/usr/bin/env python
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
# original version: 2015-06-10

import sys
from os import path
from optparse import OptionParser,OptionGroup
from pyminc.volumes.factory import *
import numpy as np
import numpy.ma as ma    # masked arrays
from scipy import ndimage,stats
import pandas as pd      # csv files
import pprint


def checkFiles(df, filenamecolumn='filename', verbose=False):
    """check files in filename column of pandas data frame df exist"""
    if filenamecolumn not in df.columns:
        raise IndexError('could not find filename column', filenamecolumn)
    filenames = list(df[filenamecolumn])
    for f in filenames:
        if not path.exists(f):
            sys.exit('Could not open file: %s' % f)
    return filenames


def getDataFrame(filename, verbose=False):
    df = pd.read_csv(filename)
    return df


def readImages(df, maskarray, filenamecolumn='filename', verbose=False):
    """reads image volumes according to filename column in pandas data frame df
    and returns a concatenated 4D volume"""
    filenames = checkFiles(df, filenamecolumn, verbose)
    vol0 = volumeFromFile(filenames[0])
    volshape = vol0.data.shape
    nvoxels = np.sum(maskarray)
    vols = np.zeros((len(filenames),nvoxels))
    for i,f in enumerate(filenames):
        vols[i,:] = volumeFromFile(f).data[maskarray].ravel()
    return vols


def createLabelVol(ar, outname, likefilename, verbose=False):
    """create label volume from array ar, saves it to disk and returns it"""
    outlabels = volumeLikeFile(likefilename, outname, dtype='ushort',
                                            volumeType="ushort", labels=True)
    outlabels.data = ar
    if verbose:
        print 'saving cluster labels to', outname
    outlabels.writeFile()
    #outlabels.closeVolume()
    return outlabels


def labeled_volume(filename, threshold, outname, maskarray, verbose=False):
    """ threshold volume from filename
    threshold is (a,b), such that <= a and >= b are kept if within maskarray"""
    if len(threshold) < 2 or threshold[0] > threshold[1]:
        raise ValueError('Threshold is not an ordered tuple')
    if verbose:
        print "thresholding with", threshold
    invol = volumeFromFile(filename)
    # threshold  (set outside threshold(s) and mask to zero)
    threshvol = invol.data.copy()
    threshvol[((invol.data > threshold[0])
              & (invol.data < threshold[1]))
              | (~maskarray) ] = 0.0
    # connected components labelling
    cluster_array, num_features = ndimage.label(threshvol)
    if verbose:
        print "found", num_features, "clusters"
    clustervol = createLabelVol(cluster_array, outname, likefilename=filename, verbose=verbose)
    return clustervol


def getCoG(labelvol, labelvalue, maskarray=None, verbose=False):
    """returns dictionary with centre of gravity in world and voxel coordinates"""
    labeldata = np.array(labelvol.data)==labelvalue
    if not isinstance(maskarray, type(None)):
        labeldata[~maskarray] = False
    # cast as np array since centre of mass doesn't like memory mapped arrays
    cog = ndimage.measurements.center_of_mass(labeldata)
    return cog


def vox2world(coord, likevol):
    """transform voxel to world coordinates"""
    return likevol.getStarts() + np.array(coord)*likevol.getSeparations()


def addCoordColumns(df, coord, columnname, dimnames, sep='-'):
    """add a column for each coordinate with suffix of likevol's dimnames"""
    for i,dim in enumerate(dimnames):
        c = columnname + sep + dim
        df[c] = coord[i]


def getValues(vols, labelvol, labelvalue, ids, maskarray=None, verbose=False):
    """Extract series of values from vols where labelvol equals labelvalue."""
    # requires vols to be pre-masked with maskarray!!!
    m = labelvol.data[maskarray]==labelvalue
    meanvalues = np.sum(vols*m, axis=1) / np.sum(m)
    df = pd.DataFrame({
         'id'    : ids
        ,'value' : meanvalues
    })
    return df


def getPeakValues(vols, peakvolname, labelvol, labelvalue, ids,
                  maskarray=None, smoothing_sigma=None, verbose=False):
    """Extract series of values from vols where labelvol equals labelvalue.
    The panda data frame df provides the ID's.
    peak values are extracted based on peaks in peakvolume """

    labelmask = (labelvol.data==labelvalue) & maskarray
    peakvol = volumeFromFile(peakvolname)

    if smoothing_sigma:
        peakvol.data = ndimage.filters.gaussian_filter(peakvol.data, smoothing_sigma)
        if verbose:
            print "peak estimation: smoothing with sigma of", smoothing_sigma

    # returns index of first extreme (regions where mask is True get masked/removed!)
    maxlinidx = ma.array(peakvol.data, mask=~labelmask).argmax()
    minlinidx = ma.array(peakvol.data, mask=~labelmask).argmin()

    fullidx = np.arange(peakvol.data.size).reshape(peakvol.data.shape)[maskarray]
    maxmaskidx = np.where(fullidx==maxlinidx)[0][0]
    minmaskidx = np.where(fullidx==minlinidx)[0][0]

    # returns first extreme (regions where mask is True get masked/removed!)
    peakidx = {
          'max' : np.unravel_index(maxlinidx, peakvol.data.shape)
        , 'min' : np.unravel_index(minlinidx, peakvol.data.shape)
    }

    peakvals = {
          'max' : vols[:,maxmaskidx]
        , 'min' : vols[:,minmaskidx]
    }

    df = pd.DataFrame()
    for extreme in ['max', 'min']:
        df2 = pd.DataFrame({
                'id' : ids
            , 'value' : peakvals[extreme]
            , 'peakvalue' : peakvol.data[peakidx[extreme]]
            , 'peaktype' : extreme
        })
        addCoordColumns(df2,           peakidx[extreme] , 'cog_vox',
                        dimnames=labelvol.dimnames)
        addCoordColumns(df2, vox2world(peakidx[extreme], likevol=peakvol), 'cog_mm' ,
                        dimnames=labelvol.dimnames)
        df = pd.concat([df,  df2], ignore_index=True)

    df['nvoxels'] = 1.0
    df['volume']  = 1.0*np.prod(peakvol.getSeparations())
    df['type']  = 'peak'

    return df


def clusterstats(vols, labelvol, df, maskarray, atlasname, defsname, voxelthresh,
                                    idcolumn, statsvolname=None, verbose=False):
    """create statistics based on df
    vols and labelvol need to be in the same space"""
    if atlasname and defsname:
        atlasvol = volumeFromFile(atlasname)
        defs = getDataFrame(defsname)

    voxelvolume = np.prod(labelvol.getSeparations())
    dfstats = pd.DataFrame()

    for i in range(1,np.max(labelvol.data)+1):
        nvoxels_cluster = np.sum(labelvol.data==i)
        if voxelthresh and nvoxels_cluster < voxelthresh:
            continue
        if verbose:
            print 'processing cluster index', i

        # get values averages across cluster
        dftmp = getValues(vols, labelvol, i, df[idcolumn].values, maskarray, verbose=verbose)
        dftmp['nvoxels'] = nvoxels_cluster
        dftmp['volume']  = nvoxels_cluster*voxelvolume
        dftmp['type']  = 'cluster'

        # add CoG
        cog = getCoG(labelvol, i, verbose=verbose)
        addCoordColumns(dftmp, cog, 'cog_vox', dimnames=labelvol.dimnames)
        addCoordColumns(dftmp, vox2world(cog, likevol=labelvol), 'cog_mm',
                        dimnames=labelvol.dimnames)

        # for enclosing structure
        if atlasname and defsname:
            # enclosing structure based on CoG
            #structure_name, labelvalue = coord_to_structure(cog, atlasvol, defs, verbose=verbose)
            # enclosing structure based on maximum overlap with cluster
            labelvalue, structure_name, side = maximum_overlap_structure((labelvol.data==i), atlasvol, defs, verbose=verbose)

            nvoxels_label = np.sum(maskarray & (atlasvol.data==labelvalue))
            # get values averages across the enclosing label
            dflabel = getValues(vols, atlasvol, labelvalue, df[idcolumn].values, maskarray, verbose=verbose)
            dflabel['nvoxels']   = nvoxels_label
            dflabel['volume']    = nvoxels_label*voxelvolume
            dflabel['type']      = 'atlaslabel'

            # add CoG
            cog = getCoG(atlasvol, labelvalue, maskarray, verbose=verbose)
            addCoordColumns(dflabel, cog, 'cog_vox', dimnames=labelvol.dimnames)
            addCoordColumns(dflabel, vox2world(cog, likevol=labelvol), 'cog_mm',
                            dimnames=labelvol.dimnames)

            dftmp = pd.concat([dftmp, dflabel])

        # for peak values
        if statsvolname:
            dfpeak = getPeakValues(vols, statsvolname, labelvol, i, df[idcolumn].values, maskarray=maskarray, smoothing_sigma=None, verbose=verbose)
            dftmp = pd.concat([dftmp, dfpeak])

        if atlasname and defsname:
            dftmp['structure'] = structure_name
            dftmp['side'] = side
            dftmp['label']     = labelvalue

        dftmp['cluster_index'] = i
        dfstats = pd.concat([dfstats, dftmp])

    return dfstats


def writeClusterstats(df, outname, selectcols=None, verbose=False):
    """writes statistics to file"""
    df.to_csv(outname, index=False)


def maximum_overlap_structure(maskarray, labelvol, defsdf, verbose=False):
    """ ... """
    labelvalue = np.bincount(labelvol.data[maskarray].astype('int')).argmax()
    tmp_right = defsdf.loc[(defsdf['right label']==labelvalue),'Structure'].values
    tmp_left  = defsdf.loc[(defsdf['left label']==labelvalue), 'Structure'].values

    if tmp_right and tmp_left:
        if tmp_right[0] != tmp_left[0]:
            raise ValueError('label', labelvalue)
        structure_name = tmp_right[0]
        side = 'both'
    elif tmp_right:
        structure_name = tmp_right[0]
        side = 'right'
    elif tmp_left:
        structure_name = tmp_left[0]
        side = 'left'
    else:
        structure_name = 'label' + str(labelvalue)
        side = 'NA'
    return (labelvalue, structure_name, side)


def coord_to_structure(coord, labelvol, defsdf, verbose=False):
    """retrieve structure label of coord"""
    #r = range(len(labelvol.dimnames))
    #idx = sorted(r,key=sorted(r,key=labelvol.dimnames.__getitem__).__getitem__)
    #coord = [coord[e] for e in [2,0,1]]
    coord = tuple([int(round(c)) for c in coord])
    labelvalue = labelvol.data[coord]     # difference between label==0 and ==None !!!
    labelvalue = int(labelvalue)
    # .values accesses underlying numpy array
    structure_name = defsdf.loc[(defsdf['right label']==labelvalue) | (defsdf['left label']==labelvalue),'Structure'].values[0]
    if not len(structure_name):
        structure_name = 'label' + str(value)
    return (structure_name, labelvalue)


def printPrettyPandas(df):
    """pretty print using whole width of terminal"""
    pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])
    pd.set_option('display.precision', 3)
    print df



###############################################################################
# MAIN
###############################################################################
if __name__ == "__main__":
    usage = """usage: %prog [-h/--help] [options] MINCFILE.."""
    description = """median filter with different kernel sizes. """
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("-m", "--mask", dest="mask",
                      help="mask name, find maximum within masked area",
                      type='string', default='')
    parser.add_option("-f", "--filetable", dest="filetable",
                      help="csv file with filename associations",
                      type='string', default="")
    parser.add_option("-i", "--inputvol", dest="inputvol",
                      help="filename of input volume that is used for thresholding",
                      type='string', default="")
    parser.add_option("-t","--absthresh",
                      help='threshold absolute input',
                      type="float", default=None)
    parser.add_option("-l","--lowerthresh",
                      help='threshold input with this lower threshold (zero below)',
                      type="float", default=None)
    parser.add_option("-u","--upperthresh",
                      help='threshold input with this upper threshold',
                      type="float", default=None)
    parser.add_option("-a", "--atlas", dest="atlas",
                      help="label volume",
                      type='string', default=None)
    parser.add_option("-d", "--defs", dest="defs",
                      help="csv file with label value - name association",
                      type='string', default=None)
    # better call this cluster size threshold
    parser.add_option("--voxelthresh", dest="voxelthresh",
                      help="only report clusters with more voxels than threshold",
                      type='float', default=None)
    parser.add_option("-o", "--outbase", dest="outbase",
                      help="output basename",
                      type='string', default=None)
    parser.add_option("--idcolumn", dest="idcolumn",
                      help="id column in table csv file",
                      type='string', default=None)
    parser.add_option("--head", dest="head",
                      help="show first n rows of output csv",
                      action="store_true", default=False)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="more verbose output",
                      action="store_true", default=False)

    (options, args) = parser.parse_args()

    if options.inputvol and not path.exists(options.inputvol):
        raise IOError('Could not open input volume: %s' % options.inputvol)

    if options.mask and path.exists(options.mask):
        maskarray = volumeFromFile(options.mask).data>0.5
    else:
        #raise IOError('Could not open mask: %s' % options.mask)
        if options.verbose:
            print 'creating whole volume mask'
        invol = volumeFromFile(options.inputvol)
        maskarray = np.ones(invol.data.shape)>0.5
        invol.closeVolume()

    if not options.filetable:
        parser.error("Filename table file not supplied")

    filedf = getDataFrame(options.filetable, verbose=options.verbose)
    vols = readImages(filedf, maskarray=maskarray, filenamecolumn='filename',
                                                        verbose=options.verbose)

    if not any((options.absthresh, options.lowerthresh, options.upperthresh)):
        threshold = (0.0, 0.0)
    elif options.absthresh:
        threshold = (-abs(options.absthresh), abs(options.absthresh))
    elif options.lowerthresh:
        threshold = (-float('inf'),options.lowerthresh)
    elif options.upperthresh:
        threshold = (options.upperthresh,float('inf'))

    if options.outbase:
        outname = options.outbase + '_index.mnc'
    else:
        outname = None

    labelvol = labeled_volume(options.inputvol, threshold,
                              outname=outname,
                              maskarray=maskarray,
                              verbose=options.verbose)

    stats = clusterstats(vols, labelvol, filedf,
                         maskarray=maskarray,
                         atlasname=options.atlas,
                         defsname=options.defs,
                         voxelthresh=options.voxelthresh,
                         idcolumn=options.idcolumn,
                         statsvolname=options.inputvol,
                         verbose=options.verbose)

    if options.outbase:
        outname = options.outbase + '.csv'
        writeClusterstats(stats, outname, selectcols=None, verbose=options.verbose)

    if options.head:
        printPrettyPandas(stats.head(30))

    labelvol.closeVolume()
    #maskvol.closeVolume()

    #dfcog = pd.DataFrame(np.array(dict(cog['voxelcoords']).values()).reshape((1,3)), columns=dict(cog['voxelcoords']).keys())
    #dfcluster.copy(deep=True)


def testing():
    """
clusterstat.py -v -f ../global/all_mice_local-fwhm0.1_short.csv -i all_anova_Cage_initgene_fdr.mnc --mask ../reg_all/hr_masks/mask.mnc -u 0.1 --voxelthresh 500 -o foo --idcolumn id --atlas atlas.mnc --defs Dorr_2008_mapping_of_labels.csv

library(RMINC)

defs <- 'tmp/tmp_defs.csv'
atlas <- 'foo_index.mnc'
table <- '../global/all_mice_local-fwhm0.1_short.csv'
t <- read.csv(table)
t$filename <- as.character(t$filename)

o <- anatGetAll(t$filename, atlas=atlas, method="means", defs=defs, side='left')
oo <- read.csv('foo.csv')

mincGetWorldVoxel(t$filename, c(-1.96, 3.18, -1.51))
mincGetVoxel(t$filename, c(48,203,77))


o[,'cluster2']
subset(oo, ROI_type=='cluster' & cluster_index==2)$value


> o[,'cluster2']
[1] -0.01421306  0.12470805 -0.03438621  0.01336526  0.05975349  0.02671909 -0.03404596 -0.12932670  0.06470247
> subset(oo, ROI_type=='cluster' & cluster_index==2)$value
[1] -0.01421306  0.12470805 -0.03438621  0.01336526  0.05975349  0.02671909 -0.03404596 -0.12932670  0.06470247


# df['x'],df['y'],df['z'] = (1,2,3)


#cog = (0,0,0)
#cog_mm = labelvol.getStarts() + np.array(cog)*labelvol.getSeparations()
#dims = labelvol.getDimensionNames()
#tmp = {'voxelcoords':zip(dims,cog),
#       'worldcoords':zip([e+unitstring for e in dims],cog_mm)}

#else:
#    ids = [str(i+1) for i in range(len(df))]
#    idcolumn = 'id'

#dftmp['max_val'] = peakvol.data[maxvoxidx]
#dftmp['min_val'] = peakvol.data[minvoxidx]
#maxidx = peakvol.getStarts() + np.array(maxvoxidx)*peakvol.getSeparations()
#minidx = peakvol.getStarts() + np.array(minvoxidx)*peakvol.getSeparations()
#addCoordColumns(dftmp, maxidx, 'max_mm', dimnames=labelvol.dimnames)
#addCoordColumns(dftmp, minidx, 'min_mm', dimnames=labelvol.dimnames)
#addCoordColumns(dftmp, maxvoxidx, 'max_vox', dimnames=labelvol.dimnames)
#addCoordColumns(dftmp, minvoxidx, 'min_vox', dimnames=labelvol.dimnames)

#if verbose: print 'peak voxel maximum', maxvoxidx, maxidx, peakvol.data[maxvoxidx]
#if verbose: print 'peak voxel minimum', minvoxidx, minidx, peakvol.data[minvoxidx]

#dfstats.rename(columns={'value': 'cluster'}, inplace=True)
    """


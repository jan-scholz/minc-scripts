#!/usr/bin/env python
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
# original version: 2015-06-10

from optparse import OptionParser,OptionGroup
from pyminc.volumes.factory import *
from scipy import ndimage,stats
import numpy as np
from os import path
import sys
import numpy.ma as ma    # masked arrays
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
    """saves label volume to disk and returns it"""
    outlabels = volumeLikeFile(likefilename, outname, dtype='ushort',
                                            volumeType="ushort", labels=True)
    outlabels.data = ar
    if verbose:
        print 'saving cluster labels to', outname
    outlabels.writeFile()
    #outlabels.closeVolume()
    return outlabels


def labeled_volume(filename, threshold, outname, maskarray, verbose=False):
    """ .. """
    if len(threshold) < 2 or threshold[0] > threshold[1]:
        raise ValueError('Threshold is not an ordered tuple')
    if verbose:
        print "thresholding with", threshold
    invol = volumeFromFile(filename)

    # threshold  (set outside threshold(s) and mask to zero)
    threshvol = invol.data.copy()
    threshvol[((invol.data >= threshold[0])
              & (invol.data <= threshold[1]))
              | (~maskarray) ] = 0.0
    # connected components labelling
    cluster_array, num_features = ndimage.label(threshvol)
    if verbose:
        print "found", num_features, "clusters"
    clustervol = createLabelVol(cluster_array, outname, likefilename=filename,
                                                                verbose=verbose)
    return clustervol


def getCoG(labelvol, labelvalue, unitstring=' (mm)', verbose=False):
    """returns dictionary with centre of gravity in world and voxel coordinates"""
    cog = ndimage.measurements.center_of_mass(labelvol.data==labelvalue)
    #cog = (0,0,0)
    cog_mm = labelvol.getStarts() + np.array(cog)*labelvol.getSeparations()
    dims = labelvol.getDimensionNames()
    tmp = {'voxelcoords':zip(dims,cog),
           'worldcoords':zip([e+unitstring for e in dims],cog_mm)}
    return tmp


def getValues(labelvol, labelvalue, vols, df, idcolumn, maskarray, verbose=False):
    """get values from vols, pandas data frame df"""
    m = labelvol.data[maskarray]==labelvalue
    bar = np.sum(vols*m, axis=1) / np.sum(m)

    if idcolumn:
        ids = df[idcolumn].values
    else:
        ids = [str(i+1) for i in range(len(df))]
        idcolumn = 'id'

    dftmp = pd.DataFrame({
         idcolumn : ids
        ,'value' : bar
    })
    return dftmp


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
        dfcluster = getValues(labelvol, i, vols, df, idcolumn, maskarray,
                                                                verbose=verbose)
        dfcluster.rename(columns={'value': 'cluster'}, inplace=True)

        dfcluster = pd.concat([dfcluster,
            pd.DataFrame({
                 'cluster_index'   : [i]
                ,'nvoxels_cluster' : [nvoxels_cluster]
                ,'volume_cluster'  : [nvoxels_cluster*voxelvolume]
                #,'type'            : ['cluster']
            }, index=dfcluster.index)], axis=1)

        cog = getCoG(labelvol, i, verbose=verbose)
        for k,v in dict(cog['voxelcoords']).iteritems():
            # replicates single values row-wise
            dfcluster[k] = pd.Series([v], index=dfcluster.index)
        for k,v in dict(cog['worldcoords']).iteritems():
            dfcluster[k] = pd.Series([v], index=dfcluster.index)

        # for enclosing structure
        if atlasname and defsname:
            # zip combines all first/second tuple elements into seperate tuples
            structure_name, labelvalue = coord_to_structure(zip(*cog['voxelcoords'])[1],
                           atlasvol, defs, verbose=verbose)
            nvoxels_label = np.sum(atlasvol.data[(maskarray) &
                                                 (atlasvol.data==labelvalue)])
            # get values averages across the enclosing label
            dflabel = getValues(atlasvol, labelvalue, vols, df, idcolumn,
                                                     maskarray, verbose=verbose)
            dflabel.rename(columns={'value': 'label'}, inplace=True)

            dfcluster = pd.concat([dfcluster,
                pd.DataFrame({
                    'structure'        : [structure_name]
                    ,'nvoxels_label' : [nvoxels_label]
                    ,'volume_label'  : [nvoxels_label*voxelvolume]
                    ,'label_value'     : [labelvalue]
                }, index=dfcluster.index)], axis=1)

            dfcluster = pd.merge(dfcluster, dflabel)
            id_vars = list(set(dfcluster.columns)-set(['cluster','label']))
            dfcluster = pd.melt(dfcluster, id_vars=id_vars, var_name='ROI_type')   #, value_name='myValname')

        dfstats = pd.concat([dfstats, dfcluster])

    return dfstats


def writeClusterstats(df, outname, selectcols=None, verbose=False):
    """writes statistics to file"""
    df.to_csv(outname, index=False)


def coord_to_structure(coord, labelvol, defsdf, verbose=False):
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
    """docstring for printPrettyPandas"""
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
    parser.add_option("--voxelthresh", dest="voxelthresh",
                      help="only report clusters with more voxels than threshold",
                      type='float', default=None)
    parser.add_option("-o", "--outbase", dest="outbase",
                      help="output basename",
                      type='string', default=None)
    parser.add_option("--idcolumn", dest="idcolumn",
                      help="id column in table csv file",
                      type='string', default=None)
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

    #printPrettyPandas(stats)

    labelvol.closeVolume()
    #maskvol.closeVolume()

    #dfcog = pd.DataFrame(np.array(dict(cog['voxelcoords']).values()).reshape((1,3)), columns=dict(cog['voxelcoords']).keys())
    #dfcluster.copy(deep=True)


def testing():
    """
clusterstat.py -v -f ../global/all_mice_local-fwhm0.1_short.csv -i all_anova_Cage_initgene_fdr.mnc --mask ../reg_all/hr_masks/mask.mnc -u 0.1 --voxelthresh 500 -o foo --idcolumn id --atlas atlas.mnc --defs Dorr_2008_mapping_of_labels.csv

library(RMINC)

defs <- 'tmp_defs.csv'
atlas <- 'foo_index.mnc'
table <- '../global/all_mice_local-fwhm0.1_short.csv'
t <- read.csv(table)
t$filename <- as.character(t$filename)

o <- anatGetAll(t$filename, atlas=atlas, method="means", defs=defs, side='left')
oo <- read.csv('foo.csv')



o[,'cluster2']
subset(oo, ROI_type=='cluster' & cluster_index==2)$value


> o[,'cluster2']
[1] -0.01421306  0.12470805 -0.03438621  0.01336526  0.05975349  0.02671909 -0.03404596 -0.12932670  0.06470247
> subset(oo, ROI_type=='cluster' & cluster_index==2)$value
[1] -0.01421306  0.12470805 -0.03438621  0.01336526  0.05975349  0.02671909 -0.03404596 -0.12932670  0.06470247
    """


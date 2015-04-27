#!/usr/bin/env python
# Calculates TFCE scores and optionally estimates quantiles of maxima
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
#         (largly based on version by Jason Lerch)
#
# TODO:
#   - restric processing to bounding box of mask?
#   - throw error when dryrun and -q option at the same time


from optparse import OptionParser,OptionGroup
from pyminc.volumes.factory import *
from scipy import ndimage,stats
from numpy import array,arange,prod
from os import remove,path


def tfce(invol, outvol, dh=0.1, E=0.5, H=2.0, mask='', verbose=False):
    outvol.data #make sure data has been loaded
    voxel_volume = prod(invol.separations)
    # step through data with dh increments
    for h in arange(0, invol.data.max(), dh):
        # threshold the data with current height
        thresh = array( invol.data > h, "uint8")
        # connected components labelling
        labeled_array, num_features = ndimage.label(thresh)
        #if verbose: print "height %d: %d features" % (h, num_features)
        # compute TFCE (array)
        if num_features > 0:
            # compute the size of each label and modulate label size by voxel volume
            sizes = array(ndimage.sum(thresh, labeled_array, range(num_features+1))) * voxel_volume
            # create an array where each voxel's value is the size of its patch
            size_array = sizes[labeled_array]
            outvol.data = outvol.data + dh*h**H * size_array**E
        else:
            # stop early if no more features
            break


def tfce_both(invol, outvol, dh=0.1, E=0.5, H=2.0, mask='', verbose=False):
    # calculate tfce scores for positive and negative values
    tfce(invol, outvol, dh, E, H, mask, verbose)
    invol.data = invol.data * -1
    tfce(invol, outvol, dh, E, H, mask, verbose)


def tfce_max(vol, mask=''):
    # get the maximum within the volume or optional mask
    if mask:
        #if verbose: print 'using mask', mask
        maskvol = volumeFromFile(mask)
        labels = array(maskvol.data > 0.5, "uint8")
    else:
        labels = None
    maxtfce = ndimage.measurements.maximum(vol.data, labels)
    return maxtfce


def run_tfce(input_filenames, suffix='_tfce',
                dh=0.1, extent=0.5, height=2.0,
                mask='',
                dryrun=False,
                no_image_output=False,
                output_max=False,
                calc_quantile=False,
                clobber=True,
                verbose=False):
    maxima = []
    for f in input_filenames:
        if verbose: print "processing %s" % f
        (root, ext) = path.splitext(f)
        output_basename = root + suffix
        max_output_filename = output_basename + '_max.txt'
        if dryrun:
            command = __file__
            for a in ['suffix', 'dh', 'extent', 'height','mask']:
                command = command + ' --' + a + '=' + str(locals()[a])
            # add --no-clobber if clobber==True
            command = command + ' --no-image-output ' + ' --save-maximum ' + f
            print command
        else:
            if calc_quantile:
                with open(max_output_filename, 'r') as f:
                    tmp =  float(f.readline())
            else:
                if clobber or not path.exists(max_output_filename):
                    invol = volumeFromFile(f)
                    outvol = volumeFromInstance(invol, output_basename + '.mnc', dtype='ushort')
                    tfce_both(invol, outvol, dh, extent, height, mask, verbose)
                    tmp = tfce_max(outvol, mask=mask)
                    if not no_image_output:
                        outvol.writeFile()
                    else:
                        remove(output_basename + '.mnc')         # dangerous workaround until volumeFromInstance doesn't create corrupt files on disk anymore
                    outvol.closeVolume()
                    if output_max:
                        f = open(max_output_filename, 'w')
                        f.write('%0.8f\n' % tmp)
                        f.close()
                else:
                    if verbose: print "using %s, unset '--no-clobber' option to overwrite" % max_output_filename
                    with open(max_output_filename, 'r') as f:
                        tmp =  float(f.readline())
            maxima.append(tmp)
    if len(maxima)>1:
        max_tfce_quantiles(maxima)


def max_tfce_quantiles(values, probs=[0.9,0.95,0.99]):
    print '# quantiles from %d maxima' % len(values)
    print 'threshold', 'quantile'                               # matches R type 7 interpolation
    for p,q in zip(probs, stats.mstats.mquantiles(values, probs, alphap=1, betap=1)):
        print '%9.2f %8.2f' % (p, q)


###############################################################################
# MAIN
###############################################################################
if __name__ == "__main__":
    usage = """usage: %prog [-h/--help] [options] MINCFILE.."""
    description = """Applies the threshold free cluster enhancement (TFCE) to a statistics image. """
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("-d", "--dh", dest="dh",
                      help="Increments over which to compute TFCE [default: %default]",
                      type="float", default=0.1)
    parser.add_option("-E", "--extent", dest="E",
                      help="Power by which to raise the extent [default: %default]",
                      type="float", default=0.5)
    parser.add_option("-H", "--height", dest="H",
                      help="Power by which to raise the height [default: %default]",
                      type="float", default=2.0)
    parser.add_option("-n", "--dryrun", dest="dryrun",
                      help="output commands for parallel processing, do not run",
                      action="store_true", default=False)
    parser.add_option("-x", "--save-maximum", dest="output_max",
                      help="save save-maximum, implied when --dryrun",
                      action="store_true", default=False)
    parser.add_option("--no-image-output", dest="no_image_output",
                      help="refrain from saving TFCE images, implied when --dryrun",
                      action="store_true", default=False)
    parser.add_option("-m", "--mask", dest="mask",
                      help="mask name, find maximum within masked area",
                      type='string', default="")
    parser.add_option("-q", "--quantile", dest="calc_quantile",
                      help="calculate quantiles from previously saved text files\nuse after --dryrun commands have run",
                      action="store_true", default=False)
    parser.add_option("--suffix", dest="suffix",
                      help="suffix appended to input files for output file names",
                      type='string', default='_tfce')
    parser.add_option("--no-clobber", dest="clobber",
                      help="more verbose output",
                      action="store_false", default=True)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="more verbose output",
                      action="store_true", default=False)

    group = OptionGroup(parser, "EXAMPLE 1 (serial, may take a long time)", path.basename(__file__) + " -m mask.mnc input*.mnc")
    parser.add_option_group(group)
    group = OptionGroup(parser, "EXAMPLE 2a (parallel, pre-processing)", path.basename(__file__) + " --dryrun -m mask.mnc input*.mnc > commands; cat commands | while read c; do sge_batch $c; done")
    parser.add_option_group(group)
    group = OptionGroup(parser, "EXAMPLE 2b (aggregation of results from 2a)", path.basename(__file__) + " --quantile -m mask.mnc input*.mnc")
    parser.add_option_group(group)
    group = OptionGroup(parser, "REFERENCE", "Smith and Nichols. Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. NeuroImage (2009) vol. 44 (1) pp. 83-98")
    #group.add_option("-g", action="store_true", help="Group option.")
    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        parser.error("Incorrect number of arguments")
    
    if options.mask and not path.exists(options.mask):
        raise IOError('Could not open mask: %s' % options.mask)

    for f in args:
        if not path.exists(f):
            raise IOError('Could not open file: %s' % f)

    run_tfce(
        args,
        suffix=options.suffix,
        dh=options.dh,
        extent=options.E,
        height=options.H,
        mask=options.mask,
        dryrun=options.dryrun,
        no_image_output=options.no_image_output,
        output_max=options.output_max,
        calc_quantile=options.calc_quantile,
        clobber=options.clobber,
        verbose=options.verbose
    )


# testing
#
# TFCE      
# time for f in out0000*_mean.mnc; do echo TFCE.py $f `basename $f .mnc`_tfce1.mnc; done | parallel -j1
#  real	0m47.952s
#  user	0m37.500s
#  sys	0m7.030s

#  real	0m45.625s
#
# time for f in out0000*_tfce1.mnc; do mincstats -quiet -max_buffer_size_in_kb 10066240 -mask mask.mnc  -mask_binvalue 1 -max $f; done > TFCE1_max.txt;
#
# R --no-save --no-restore -e 'print(quantile(read.table("TFCE1_max.txt")$V1,c(0.9,0.95,0.99)))'


# TFCE2
# time TFCE2b.py -v -m mask.mnc out0000*_mean.mnc
#real	0m55.834s
#user	0m46.660s
#sys	0m4.860s

# time TFCE2b.py --dryrun -m mask.mnc out0000*_mean.mnc | parallel -v
# real	0m24.101s
# real	0m22.651s
# user	2m0.800s
# sys	0m17.740s
# time TFCE2b.py --quantile -m mask.mnc out0000*_mean.mnc
#real	0m0.965s
#user	0m0.260s
#sys	0m0.170s


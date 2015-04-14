#!/usr/bin/env python
#
# 2015-04-10
# concate two tag files in column direction for input in tagtoxfm


from optparse import OptionParser
import numpy as np
from os import path
from StringIO import StringIO


def readTagFile(tagfilename, verbose=False):
    """read tag file and return numpy array with tags
    (6 columns, a point per row)
    returns tags as numpy array"""

    with open(tagfilename, 'r') as f:
        filetype = f.readline().strip()
        nvolumes = f.readline().split('=')[1].strip().strip(';')
        lines = f.readlines()
        for i,l in enumerate(lines):
            if 'points' in l.lower():
                break
        lines = ''.join(lines[(i+1):]).rstrip().rstrip(';')

    if filetype.lower() != 'mni tag point file':
        raise ValueError('not a "MNI Tag Point File", claims to be', filetype)
    if nvolumes != '1':
        raise ValueError('reading of only 1 volume supported, not:', nvolumes)

    f = StringIO(lines)
    tags = np.loadtxt(f)

    if verbose:
        print 'read %i points from tag file %s' % (tags.shape[0],tagfilename)

    return tags


def cbind_tagfiles(outname, filename1, filename2, verbose=False):
    """read two tagfiles and column-bind them
    (keeps only coordinates, last 3 coordinates get discarded)
    returns combined tags numpy array"""

    tags1 = readTagFile(filename1, verbose)
    tags2 = readTagFile(filename2, verbose)
    if tags1.shape[0] != tags2.shape[0]:
        raise ValueError('input files have different number of rows')
    outtags = np.hstack((tags1[:,0:3],tags2[:,0:3]))

    footer = ';'
    header = """MNI Tag Point File
Volumes = 2;
%% File created automatically by %s. Input tag files:
%% %s, %s

Points = """ % (path.basename(__file__), filename1, filename2)

    np.savetxt(outname, outtags, fmt='%g', header=header, footer=footer, comments='')

    return outtags


def main():
    usage = "usage: %prog -o OUTPUT TAGFILE1 TAGFILE2"
    description = "smooth labels through prob. exchange"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--output", "-o", dest="output",
                      help="output file",
                      type="string")
    parser.add_option("-v", "--verbose",
                      help='verbose output', action='store_true', default=False)
    (options, args) = parser.parse_args()

    if (len(args) < 2):
        parser.error("require exactly two tag files as input")

    if options.output is None:
        parser.error("output is missing")

    for f in args:
        if not path.exists(f):
            raise IOError('could not find file: ' + f)

    cbind_tagfiles(options.output, args[0], args[1], verbose=options.verbose)


if __name__ == "__main__":
    main()


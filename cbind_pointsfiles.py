#!/usr/bin/env python
#
# 2015-04-10
# concate two tag or obj files in column direction for input in tagtoxfm


from optparse import OptionParser
import numpy as np
from os import path
from StringIO import StringIO


def readObjFile(filename, verbose=False):
    """read obj file and return numpy array with points
    returns numpy array"""
    with open(filename, 'r') as f:
        # header line
        headerline = f.readline().split()
        if headerline[0].upper() != 'P':
            raise ValueError('not an OBJ file, \
                    expected first line to start with "P"', ' '.join(headerline))
        npoints = int(headerline[-1])
        lines = f.readlines()
    # extract npoints rows, the rest are normals and connectivity information
    f = StringIO(''.join(lines[0:npoints]))
    tags = np.loadtxt(f)
    if verbose:
        print 'read %i (header: %i) points from tag file %s' % \
                                            (tags.shape[0], npoints, filename)
    return tags


def readTagFile(filename, verbose=False):
    """read tag file and return numpy array with tags
    returns tags as numpy array"""
    with open(filename, 'r') as f:
        # header
        filetype = f.readline().strip()
        nvolumes = f.readline().split('=')[1].strip().strip(';')
        lines = f.readlines()
        for i,l in enumerate(lines):
            if 'points' in l.lower():
                break
        lines = ''.join(lines[(i+1):]).rstrip().rstrip(';')
    # verify header information
    if filetype.lower() != 'mni tag point file':
        raise ValueError('not a "MNI Tag Point File", claims to be', filetype)
    if nvolumes != '1':
        raise ValueError('reading of only 1 volume supported, not:', nvolumes)
    f = StringIO(lines)
    tags = np.loadtxt(f)
    if verbose:
        print 'read %i points from tag file %s' % (tags.shape[0],filename)
    return tags


def columnbind_files(outname, filenames, verbose=False):
    """read two tag/obj files and column-bind them
    (keeps only coordinates, last 3 coordinates (i.e. columns) get discarded)
    returns combined tags numpy array"""
    if len(filenames) < 2:
        raise ValueError('need at least two inputs')
    tags = []
    for f in filenames[:2]:
        informat = path.splitext(f)[1].strip('.')
        if informat=='tag':
            tags.append(readTagFile(f, verbose))
        if informat=='obj':
            tags.append(readObjFile(f, verbose))
    if tags[0].shape[0] != tags[1].shape[0]:
        raise ValueError('input files have different number of rows')
    outtags = np.hstack((tags[0][:,0:3],tags[1][:,0:3]))
    # prepare header/footer for tag file output
    footer = ';'
    header = """MNI Tag Point File
Volumes = 2;
%% File created automatically by %s. Input tag files:
%% %s, %s

Points = """ % (path.basename(__file__), filenames[0], filenames[1])
    np.savetxt(outname, outtags, fmt='%g', header=header, footer=footer, comments='')
    if verbose:
        print 'saved to file:', outname
    return outtags


def main():
    usage = "usage: %prog -o OUTPUT INFILE1 INFILE2"
    description = "column-bind TAG or OBJ files for processing with tag2xfm"
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

    columnbind_files(options.output, (args[0], args[1]), verbose=options.verbose)


if __name__ == "__main__":
    main()


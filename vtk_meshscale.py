#!/usr/bin/env python
#
# Apply vertex-wise scaling to a mesh, e.g. eigenvector from PCA analysis
#
# author: jan.scholz@mouseimaging.ca
# original version: 2015-05-15

from os import path
from optparse import OptionParser
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk


class MyParser(OptionParser):
    """alow adding usage examples in epilog"""
    def format_epilog(self, formatter):
        return '\n' + self.epilog + '\n'


def readMeshFile(filename, verbose=False):
    """Read mesh file.
    The input format is determined by file name extension."""

    informat = path.splitext(options.infilename)[1].strip('.')
    # set reader based on filename extension
    if informat=='stl':
        reader = vtk.vtkSTLReader()
    elif informat=='vtk':
        reader = vtk.vtkPolyDataReader()
    elif informat=='obj':
        reader = vtk.vtkMNIObjectReader()
    #elif informat=='tag':
    #    reader = vtk.vtkMNITagPointReader()
    else:
        raise ValueError('cannot read input format: ' + informat)
    reader.SetFileName(filename)
    reader.Update()

    if verbose:
        # nPolys*3 == R's readSTL nrow
        print "read %i polygons from file %s" % \
                               (reader.GetOutput().GetNumberOfPolys(), filename)

    return reader.GetOutput()


def readRData(filename, objectname, verbose=False):
    import rpy2.robjects as robjects
    robjects.r['load'](filename)
    a = robjects.r[objectname]
    d = dict(zip(a.names, list(a))) # only one level deep !!!
    arrayfieldnames = ['mshape', 'pcar']
    rd = dict()
    for s in arrayfieldnames:
        rd[s] = np.ascontiguousarray(d[s], dtype=np.float32)
        if verbose:
            print 'read', s, 'with dimensions', rd[s].shape
    return rd


def scaleMesh(polydata, rdict, col=1, scale=1.0, verbose=False):
    """scale vertices/points in a mesh
    pcar: columns are eigenvectors (PCs)
    col: PC number, starts at 1"""

    if col < 1:
        raise ValueError('col needs to be larger than 0')
    col = col -1; # python matrix starts at 0

    p = polydata.GetPoints()  # returns vtk points
    #d = p.GetData()          # returns vtk float array  # using mshape instead
    #npdata = vtk_to_numpy(d)                            # using mshape instead

    pcar = rdict['pcar']
    mshape = rdict['mshape']
    pcar = pcar[:,col].reshape((mshape.shape[::-1])).transpose()

    if verbose:
        print 'adding PC', col, 'with scale', scale
    mshape = mshape + scale * pcar
    vtkdata = numpy_to_vtk(num_array=mshape, deep=True, array_type=vtk.VTK_FLOAT)

    p.SetData(vtkdata)
    polydata.SetPoints(p)
    return polydata


def writeMeshFile(polydata, filename, binary=True, verbose=False):
    """Write mesh file.
    The output format is determined by file name extension. Files can be written
    in binary (default) and ASCII format."""

    outformat = path.splitext(options.outfilename)[1].strip('.')
    # set writer based on filename extension
    if outformat=='stl':
        writer = vtk.vtkSTLWriter()
    elif outformat=='vtk':
        writer = vtk.vtkPolyDataWriter()
    elif outformat=='obj':
        writer = vtk.vtkMNIObjectWriter()
    elif outformat=='tag':
        writer = vtk.vtkMNITagPointWriter()
    else:
        raise ValueError('cannot write output format' + outformat)
    #writer.SetInputConnection(polydata.GetOutputPort())
    writer.SetInputData(polydata)

    if outformat!='tag':
        if binary:
            if verbose: print 'setting ouptut to binary'
            writer.SetFileTypeToBinary()
        else:
            if verbose: print 'setting ouptut to ascii'
            writer.SetFileTypeToASCII()

    writer.SetFileName(filename)
    err = writer.Write()
    if err != 1:
        raise IOError('failed to write')

    if verbose:
        print "wrote", filename
    pass


if __name__ == "__main__":
    usage = """usage: %prog [-h/--help] -i INFILE -r VECTOR -s SCALE -o OUTFILE"""

    description = """Convert between mesh file formats.
    Currently supports reading and writing of STL, VTK, OBJ (BIC object)"""
    epilog = "Example:\n  " + \
        path.basename(__file__) + " -v -i foo.vtk -r vec.txt -s 2 -o bar.stl"

    parser = MyParser(usage=usage, description=description, epilog=epilog)

    parser.add_option("-i", "--input", dest="infilename",
                      help="no help",
                      type='string', default="")
    parser.add_option("-o", "--output", dest="outfilename",
                      help="no help",
                      type='string', default="")
    parser.add_option("--rdata", dest="rdata",
                      help="no help",
                      type='string', default="")
    parser.add_option("--objectname", dest="objectname",
                      help="no help",
                      type='string', default="")
    parser.add_option("-c", "--column", dest="column",
                      help="no help",
                      type='int', default=1)
    parser.add_option("-s", "--scale", dest="scale",
                      help="no help",
                      type='float', default=0.0)
    parser.add_option("-a", "--ascii", dest="binary",
                      help="save in ascii format",
                      action="store_false", default=True)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="more verbose output",
                      action="store_true", default=False)

    (options, args) = parser.parse_args()

    if options.infilename is '':
        parser.error('INFILE not specified (-i)')
    if options.outfilename is '':
        parser.error('OUTFILE not specified (-o)')
    if options.rdata is '':
        parser.error('RData file not specified (--rdata)')
    if options.objectname is '':
        parser.error("objectname not specified (--objectname 'LIST$FIELD)'")
    if not path.exists(options.infilename):
        parser.error('could not find INFILE')
    if not path.exists(options.rdata):
        parser.error('could not find RData file')
    outpath = path.dirname(options.outfilename)
    if outpath and not path.exists(outpath):
        parser.error('output directory does not exist: ' + outpath)


    pd = readMeshFile(options.infilename, verbose=options.verbose)
    rdict = readRData(options.rdata, options.objectname,
                                                        verbose=options.verbose)
    pd = scaleMesh(pd, rdict=rdict, col=options.column, scale=options.scale,
                                                        verbose=options.verbose)
    writeMeshFile(pd, options.outfilename, binary=options.binary,
                                                        verbose=options.verbose)




# vtk references
# https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
# http://stackoverflow.com/questions/6684306/how-can-i-read-a-vtk-file-into-a-python-datastructure


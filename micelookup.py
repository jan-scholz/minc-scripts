#!/usr/bin/env python
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
#
# TODO:
#   - dimension order
#   - world coordinates
#   - nearest label if label==0
#   - CoG of structure


#defaulttable = '/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv'
defaulttable = '/home/jscholz/resources/segmentation/c57_brain_atlas_labels.csv'

tablefiles = ['/projects/mice/share/mouse-brain-atlases/Dorr_2008_Steadman_2013_Ullmann_2013/Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv','/projects/mice/share/mouse-brain-atlases/Dorr_2008/ex-vivo/Dorr_2008_mapping_of_labels.csv','/projects/mice/share/mouse-brain-atlases/Dorr_2008/in-vivo-t2/Dorr_2008_on_in-vivo-t2_mapping_of_labels.csv','/projects/mice/share/mouse-brain-atlases/Dorr_2008/Dorr_2008_mapping_of_labels.csv']

from optparse import OptionParser
from pyminc.volumes.factory import *
import numpy as np
from os import path
import sys
import string


def getTable(tablefilename, delimiter = ",", verbose=False):
    with open(tablefilename) as f:
        tmp = f.read().splitlines()
    tmp = [e.split(delimiter) for e in tmp]
    coldiff = len(tmp[1])-len(tmp[0])         #raise ValueError('non-matching header')
    table = list()
    for i,r in enumerate(tmp[1:]):
        table.append(dict(zip(tmp[0],r[coldiff:])))
    if verbose: print 'using table %s' % (tablefilename,len(tablefilename))
    return table


# df = pd.read_csv(tablefilename)
# df.columns
# df['Structure']
# df.loc[1,'Structure']   # starts at one


def label_to_structure(value, table, verbose=False):
    if verbose: print 'label to structure'
    tmp = [e for e in table if value in (e['left label'],e['right label'])]
    if not len(tmp):
        tmp = [{'label':value}]
    return tmp

def structure_to_label(structure, table, verbose=False):
    if verbose: print 'structure to label'
    tmp = [e for e in table if string.find(e['Structure'], structure)>-1]
    #width = max([len(e['Structure']) for e in tmp])   #.ljust(width)
    if not len(tmp):
        tmp = [{'Structure':structure}]
    return tmp

def isCoord(coords):
    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False
    result = sum([is_number(e) for e in coords.split(',')])
    return result==3

def coord_to_label(coord, labelfile, inworldcoords=False, verbose=False):
    if verbose: print 'coord to label: ', coord
    vol = volumeFromFile(labelfile)
    r = range(len(vol.dimnames))
    idx = sorted(r,key=sorted(r,key=vol.dimnames.__getitem__).__getitem__)
    coord = [coord[i] for i in idx]
    if inworldcoords:
        coord = (np.array(coord)-np.array(vol.starts))/np.array(vol.separations)
    coord = tuple([int(round(c)) for c in coord])
    #print '@', coord
    label = vol.data[coord]     # difference between label==0 and ==None !!!
    #print 'l', label
    vol.closeVolume()
    return str(int(label))

def coord_to_structure(coord, labelfile, table, inworldcoords=False, verbose=False):
    if verbose: print 'coord to structure: '
    l = coord_to_label(coord, labelfile, inworldcoords=inworldcoords, verbose=verbose)
    tmp = label_to_structure(l, table, verbose=verbose)
    tmp[0]['coordinates']=coord
    return tmp

def printResults(d, onlystructure=False):
    def printNiceCoord(coord):
        return ','.join(format(f, '.2f').rstrip('0').rstrip('.') for f in coord)
    if len(d):
        for e in d:
            if len(e) < 3:
                k,v = e.iteritems().next()
                if k == 'coordinates':
                    v = '<' + printNiceCoord(v) + '>'
                print 'no results found for %s %s' % (k,v)
            else:
                s = '%s' % e['Structure']
                if not onlystructure:
                    s += ', left label: %3s, right label: %3s' % (e['left label'], e['right label'])
                    if 'coordinates' in e.keys():
                        s += ', at: <%s>' % printNiceCoord(e['coordinates'])
                print s
    else:
        sys.stderr.write('no results to print\n')


###############################################################################
# MAIN
###############################################################################
if __name__ == "__main__":
    usage = """usage: %prog [-h/--help] [options] LABEL_VALUE|STRUCTURE_NAME|LABEL_FILE|.."""
    description = """Looks up name and labels given a (partial) name, a label, or coordinates.
Use '--' to separate negative coordinates from options, e.g. %prog -l labels.mnc -- -1.2,0,0"""
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("-t", "--table", dest="table",
                      help="table file name, csv",
                      type='string', default=defaulttable)
    parser.add_option("-l", "--labelfile", dest="labelfile",
                      help="minc atlas file with labels",
                      type='string', default=None)
    parser.add_option("-w", "--worldcoords", dest="inworldcoords",
                      help="coordinates are given in world coordinates, e.g. mm",
                      action="store_true", default=False)
    parser.add_option("-s", "--only-structures", dest="onlystructure",
                      help="output only the name of the structure",
                      action="store_true", default=False)
    parser.add_option("--list", dest="listtablefiles",
                      help="list table files",
                      action="store_true", default=False)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="more verbose output",
                      action="store_true", default=False)

    (options, args) = parser.parse_args()


    if options.listtablefiles:
        print '\n'.join(tablefiles)
        sys.exit(0)

    if not len(args):
        parser.print_usage()
        sys.exit(1)

    table = getTable(options.table, verbose=options.verbose)

    for a in args:
        if a.isdigit():
            out = label_to_structure(a, table, verbose=options.verbose)
            printResults(out, onlystructure=options.onlystructure)
        elif isCoord(a):
            if not options.labelfile:
                raise TypeError('LABELFILE not specified')
            coord = [float(e) for e in a.split(',')]
            out = coord_to_structure(coord, options.labelfile, table, inworldcoords=options.inworldcoords, verbose=options.verbose)
            printResults(out, onlystructure=options.onlystructure)
        else:
            out = structure_to_label(a, table, verbose=options.verbose)
            printResults(out, onlystructure=options.onlystructure)



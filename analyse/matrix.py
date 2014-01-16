
from matplotlib import pyplot as plt

import argparse
import json
import os.path
import re


parser = argparse.ArgumentParser(description='aggregates several .json files and produces matrix comparing classification with reference')

parser.add_argument('-v', '--verbose', action='count',
        help='print some information to stderr')

parser.add_argument('-n', '--nothing', action='store_true',
        help='don\'t actually output (draw) anything')

parser.add_argument('-o', '--output', nargs=1,
        help='output file of produced graphic')

parser.add_argument('-r', '--reference', nargs=1,
        help='specify .json file that acts as a reference; dictionary indexed by filename without extension')

parser.add_argument('-f', '--filter', nargs=1, default=['.*'],
        help='filter entries to display in hitmap (python regular expression)')

parser.add_argument('-d', '--dpi', nargs=1, default=[150], type=int,
        help='specify resolution of output image; set higher resolution to get "smaller" labels')
parser.add_argument('-x', '--width', nargs=1, default=[30],
        help='set width of figure (in inches)')
parser.add_argument('-y', '--height', nargs=1, default=[20],
        help='set height of figure (in inches)')
parser.add_argument('-u', '--fractionx', nargs=1, default=[1.],
        help='fraction of image used by matrix (remainder goes to labels)')
parser.add_argument('-w', '--fractiony', nargs=1, default=[.8],
        help='fraction of image used by matrix (remainder goes to labels)')


parser.add_argument('command', choices=['lineage','hitmap', 'spoligo', 'resistance'],
        help='what kind of analysis to perform')
parser.add_argument('json', nargs='+',
        help='name of .json files to be aggregated')

args = parser.parse_args()


reference = None
if args.reference:
    f = open(args.reference[0])
    reference = json.load(f)

jsons={}
for fname in args.json:
    f = open(fname)
    key = os.path.splitext(os.path.basename(fname))[0]
    jsons[key] = json.load(f)

    assert jsons[key]['info']['format'] == 'kvarq'


def plot_matrix(M, xlabels, ylabels):
    if args.nothing:
        return
    width = args.width[0]
    height = args.height[0]
    dpi = args.dpi[0]
    plt.figure(figsize=(width, height), dpi=dpi)
    fx = args.fractionx[0]
    fy = args.fractiony[0]
    m = 0.05
    ax = plt.axes([1-fx +m,1-fy +m,fx -2*m,fy -2*m]) #TODO adjust automatically
    #print str([1-fx +m,1-fy +m,fx -2*m,fy -2*m])
    #ax = plt.axes([.05,.2,.9,.75]) #TODO adjust automatically
    #ax = plt.axes([0.05, 0.24999999999999994, 0.9, 0.7000000000000001])
    plt.imshow(matrix, interpolation='nearest', axes=ax)
    plt.xticks(range(len(xlabels)), xlabels, rotation='vertical')
    plt.yticks(range(len(ylabels)), ylabels)


def equal(lineage, reflineage):
    if lineage == reflineage:
        return True
    ls = lineage.partition('/')
    rls= reflineage.partition('/')
    if ls[0] == rls[0] and len(ls[2])*len(rls[2])==0:
        return True
    return False


def oct2bin(ostr):
    assert len(ostr)==15
    spol42 = ostr[-1]
    value = int(ostr[:-1], 8)
    if spol42=='0':
        ret = [0]
    else:
        ret = [1]
    for i in range(42):
        ret = [value%2] + ret
        value >>= 1
    return ret


if args.command == 'lineage':

    total = matches = mismatches = 0

    lineages = set()
    for name,data in jsons.items():
        #mainlineage = data['analyses']['lineage'].partition('/')[0]
        lineage = data['analyses']['lineage']
        lineages.add(lineage)
    if reference:
        for name,data in reference.items():
            if 'lineage' in data:
                lineages.add(data['lineage'])
    lineages = list(lineages)
    lineages.sort()

    def row_from_lineage(lineages, lineage):
        row = [0] * len(lineages)
        row[lineages.index(lineage)] = 1
        return row

    matrix = []
    row_labels = []
    for name,data in jsons.items():

        row_labels.append(name)
        #mainlineage = data['analyses']['lineage'].partition('/')[0]
        code = data['analyses']['spoligo']
        matrix.append(row_from_lineage(lineages, lineage))

        total += 1

        if reference and name in reference and 'lineage' in reference[name]:
            reflineage = reference[name]['lineage']
            matrix.append(row_from_lineage(lineages, reflineage))
            row_labels.append('*')

            if equal(lineage, reflineage):
                matches += 1
            else:
                mismatches += 1
                if args.verbose>1: print '%s : %s instead of %s'%(
                            name, lineage, reflineage)

    plot_matrix(matrix, lineages, row_labels)

    if args.verbose:
        try:
            print '%d strains -> %d %% matches (%d vs %d)'%(
                total, int(100.*matches/(matches+mismatches)), matches, mismatches)
        except ZeroDivisionError, e:
            print '%d strains (no refernce information)'% total


if args.command == 'spoligo':

    total = matches = mismatches = 0

    anamax = 0
    matrix = []
    row_labels = []
    for name,data in jsons.items():

        code = data['analyses']['spoligo']
        spoligo = oct2bin(code)

        analog=[0]*43
        for x,hits in data['hits'].items():
            m = re.match('^spoligo(\\d+)', x)
            if m:
                v = len(hits.split(','))
                anamax = max(anamax, v)
                analog[int(m.group(1))] = v

        total += 1

        if reference and name in reference and 'spoligo' in reference[name]:
            refcode = reference[name]['spoligo']
            refspoligo = oct2bin(refcode)

            # row 1 : matches
            matrix.append(analog)
            row_labels.append(name)
            # row 2 : spoligo
            matrix.append(spoligo)
            row_labels.append(name)
            # row 3 : reference
            matrix.append(refspoligo)
            row_labels.append('*')

            for i,x in enumerate(spoligo):
                if x == refspoligo[i]:
                    matches += 1
                else:
                    mismatches += 1

    for i in range(len(matrix)/3):
        row = matrix[i*3]
        for j in range(len(row)):
            row[j] /= float(anamax)

    if not matrix:
        print '*** no reference data tound -> matrix cannot be plotted'
    else:
        plot_matrix(matrix, ['']*43, row_labels)

    if args.verbose:
        print '%d strains -> %d %% matches (%d vs %d)'%(
            total, int(100.*matches/(matches+mismatches)), matches, mismatches)


if args.command == 'hitmap':

    filterre = re.compile(args.filter[0])
    hits = set()
    for name,data in jsons.items():
        for hit,poss in data['hits'].items():
            if filterre.match(hit):
                hits.add(hit)
    hits = list(hits)
    hits.sort()

    matrix = []
    row_labels = []
    for name,data in jsons.items():

        row = [0] * len(hits)
        for hit,poss in data['hits'].items():
            if filterre.match(hit):
                row[hits.index(hit)] = len(poss.split(','))

        matrix.append(row)
        row_labels.append(name)

    plot_matrix(matrix, hits, row_labels)


if args.output and not args.nothing:
    plt.savefig(args.output[0], dpi=args.dpi[0])
    plt.close('all')


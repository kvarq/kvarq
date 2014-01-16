
import argparse
import json
import sys
import time

parser = argparse.ArgumentParser(description='converts .json from one version to another')

parser.add_argument('-v', '--verbose', action='store_true',
        help='print some information to stderr')

parser.add_argument('json',type=file,
        help='name of .json file to process')

args = parser.parse_args()

data = json.load(args.json)


if not 'info' in data:
    # convert version 0 -> version 1

    if args.verbose: sys.stderr.write('converting '+args.json.name+'\n')

    for test,posstr in data['hits'].items():
        poss = []
        lpos = 0
        overflow = 0
        for pos in map(int, posstr.split(',')):
            if lpos<0 and pos>=0:
                sys.stderr.write('*** probable overflow in %s : %d -> %d\n'%(test,lpos,pos))
                overflow += 1
            lpos = pos
            if pos<0:
                pos = (pos-(-2147483648)) + 0x80000000
            poss.append(pos)

        data['hits'][test] = ','.join(map(str,poss))

    data['analyses'] = {
        'spoligo': data['spoligo'],
        'lineage': data['lineage'],
        }

    del data['spoligo']
    del data['lineage']

    descr = 'converted 0->1 on '+time.asctime(time.localtime())
    if overflow>0:
        descr+='; probably contains overflows that were not fixed'

    data['info'] = {
        'format':'kvarq',
        'fastq':'?',
        'when':'?',
        'descr': descr,
    }


print json.dumps(data, indent=2)


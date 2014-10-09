
import collections, csv, sys, argparse
import xlrd

parser = argparse.ArgumentParser(
        description='extract SL1, SL2, SL3 SNPs from supplementary tables')
parser.add_argument('output_format', choices=['csv', 'python'],
        help='in what format to output the data; "csv" outputs a table '
        'that shows the respective frequencies for every SNP; "python" '
        'directly outputs the Test(SNP(...)) code to be used with KvarQ')
args = parser.parse_args()

wb = xlrd.open_workbook('Table S2 samples.xlsx')
ws = wb.sheet_by_name('Samples')
sampleids = set()
for row in range(1, ws.nrows):
    sampleid = str(ws.cell(row, 1).value)
    # we're not interested in multiple samples from same individual
    sampleids.add(sampleid.split('.')[0])

wb = xlrd.open_workbook('Table S4 2014_specific_snps.xlsx')
ws = wb.sheet_by_name('2014_specific_snps')

SNP = collections.namedtuple('SNP', ('pos', 'orig', 'base'))

snps = []
for row in range(1, ws.nrows):

    pos = int(ws.cell(row, 0).value)
    orig = ws.cell(row, 1).value
    base = ws.cell(row, 2).value

    snps.append(SNP(pos, orig, base))

snps.sort(key=lambda x: x.pos)

matrix = {}
bases = 'ACGT-'

ident = pos = snpi = None
for line in file('ebov.mafft.fasta'):

    line = line.strip()

    if line[0] == '>':
        ident = line[1:]
        pos = 1
        snpi = 0
        continue

    while snpi < len(snps) and snps[snpi].pos < pos:
        snpi += 1

    while snpi < len(snps) and snps[snpi].pos < pos + len(line):

        snp = snps[snpi]
        base = line[snp.pos - pos]
        if snp.pos not in matrix:
            matrix[snp.pos] = dict([(b, []) for b in bases])
        matrix[snp.pos].setdefault(base, []).append(ident)

        snpi += 1

    pos += len(line)


if args.output_format == 'csv':
    w = csv.writer(sys.stdout)
    header = ['pos']
    for b in bases:
        header += [b + '.sierraleone', b + '.guinea', b + '.previous']

    w.writerow(header)

if args.output_format == 'python':
    print('SNPs = [')

for snp in snps:
    row = [snp.pos]
    SL = [None] * 3
    orig = None
    for b in bases:
        l = matrix[snp.pos][b]

        n_2014 = sum([1 for ident in l if '2014' in ident])
        n_sierraleone = sum([1 for ident in l if ident.split('_')[-1] in sampleids])
        n_guinea = n_2014 - n_sierraleone
        n_previous = len(l) - n_2014

        if args.output_format == 'csv':
            row += [n_sierraleone, n_guinea, n_previous]

        if args.output_format == 'python':
            if n_sierraleone == 78:
                SL[0] = b
            elif n_sierraleone == 72:
                SL[1] = b
            elif n_sierraleone == 44:
                SL[2] = b
            if n_guinea == 3 and n_previous == 20:
                orig = b

    for i in range(3):
        if SL[i] and orig and SL[i] != orig:
            print(' '*8 + "Test(SNP(genome=EBOV76, pos=%d, orig='%s', base='%s'), SL%d, gire14)," %
                    (snp.pos + 1, orig, SL[i], i+1))

    if args.output_format == 'csv':
        w.writerow(row)

if args.output_format == 'python':
    print(' '*4 + ']')


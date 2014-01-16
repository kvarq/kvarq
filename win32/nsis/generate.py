# generates file lists
# call this after building py2exe

import os, os.path, glob, shutil
from os.path import pardir

dst = os.path.join(os.path.dirname(__file__), pardir, pardir, 'dist')
inc = os.path.join(os.path.dirname(__file__), 'include')
exes = [os.path.basename(exe).lower() for exe in glob.glob(os.path.join(dst, '*.exe'))]
print exes
assert 'simple.exe' in exes

instf = open(os.path.join(inc, 'install.files'), 'w')
uninstf = open(os.path.join(inc, 'uninstall.files'), 'w')

for cand in glob.glob(os.path.join(dst, '*')):

    if cand.lower().endswith('.nsi'):
        continue

    basecand = cand[len(dst)+1:]

    if os.path.isdir(cand):
        instf.write('  File /r "%s"\n'% basecand)
        uninstf.write('  RMDir /r "$INSTDIR\\%s"\n'% basecand)
        print 'added %s\\'% basecand
    else:
        instf.write('  File "%s"\n'% basecand)
        uninstf.write('  Delete "$INSTDIR\\%s"\n'% basecand)
        print 'added %s'% basecand

instf.close()
uninstf.close()

for fname in ['installer.nsi']:
    nsi = os.path.join(os.path.dirname(__file__), fname)
    shutil.copyfile(nsi, os.path.join(dst, fname))


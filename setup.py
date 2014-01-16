from setuptools import setup, Extension
import os, os.path
import glob
import os
import sys
import shutil

import kvarq

def shoutout(msg):
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        print '\033[91m'
    print '*' * (len(msg) + 2*5)
    print '**** ' + msg + ' ****'
    print '*' * (len(msg) + 2*5)
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        print '\033[m'

kwargs = dict(options={})
additional_datafiles = []

# optional dependencies {{{1
# include (optional) xlwt, xlrd if locally installed
includes = []
try:
    import xlwt
    includes.append('xlwt')
except ImportError:
    shoutout('building app/exe file without xlwt')
try:
    import xlrd
    includes.append('xlrd')
except ImportError:
    shoutout('building app/exe file without xlrd')

# configure compilation extensions {{{1
if sys.platform == 'win32':
    # see http://stackoverflow.com/questions/2817869/error-unable-to-find-vcvarsall-bat
    if 'VS100COMNTOOLS' in os.environ: os.environ['VS90COMNTOOLS'] = os.environ['VS100COMNTOOLS']
    if 'VS110COMNTOOLS' in os.environ: os.environ['VS90COMNTOOLS'] = os.environ['VS110COMNTOOLS']
    # download pthreads from http://sourceware.org/pthreads-win32/
    pthread32path = os.path.join(os.path.dirname(__file__), 'win32', 'pthread')
    os.environ['PATH'] += os.path.pathsep + pthread32path
    engine_ext = Extension('kvarq.engine',glob.glob(os.path.join('csrc','*.c')),
            libraries=[os.path.join(pthread32path, 'pthreadVC2')], extra_compile_args=['/I', pthread32path])

else:
    engine_ext = Extension('kvarq.engine',glob.glob(os.path.join('csrc','*.c')))


# package builders {{{1

# windows : py2exe {{{2
if sys.platform == 'win32':
    kwargs['windows'] = [{
            'script': os.path.join('kvarq', 'gui', 'main.py'),
            'icon_resources': [(1, os.path.join('res', 'logo', 'TPH_DNA.iconset', 'icon_64x64.ico'))],
        }]
    kwargs['console'] = [os.path.join('kvarq', 'cli.py')]
    kwargs['py2exe'] = {
            'icon': os.path.join('res', 'logo', 'TPH_DNA.iconset', 'icon_64x64.ico'),
        }
    kwargs['options']['py2exe'] = {
            'dll_excludes': ['MSVCP90.dll'],
            'includes': includes,
        }
    additional_datafiles += [('res', [os.path.join('res', 'TPH_DNA.ico')])]

    try:
        # monkey patch py2exe command to rename executables
        import py2exe # pylint: disable=F0401
        from py2exe.build_exe import py2exe as py2exe_command # pylint: disable=F0401
        py2exe_command_run = py2exe_command.run
        def py2exe_mv_files(self):
            py2exe_command_run(self)
            print '*** renaming .exe files ***'
            for src, dst in (('cli.exe', 'kvarq.exe'), ('main.exe', 'kvarq-gui.exe')):
                src = os.path.join(self.dist_dir, src)
                dst = os.path.join(self.dist_dir, dst)
                print '%s -> %s'%(src, dst)
                if os.path.exists(dst): os.unlink(dst)
                os.rename(src, dst)
        py2exe_command.run = py2exe_mv_files

    except ImportError, e:
        pass

# OS X : py2app {{{2
if sys.platform == 'darwin':
    kwargs['options']['py2app'] = dict(
            iconfile='res/logo/TPH_DNA.icns',
            includes=includes
        )

# warning if documentation outdated {{{1
newest_rst = sorted(
        [fname for fname in glob.glob(os.path.join('docs', '*.rst'))
                if os.path.isfile(fname)],
        key=os.path.getmtime)[-1]
htmldir = os.path.join('docs', '_build', 'html')
if os.path.isdir(htmldir):
	newest_html = sorted(
		[fname for fname in glob.glob(os.path.join(htmldir, '*.html'))
			if os.path.isfile(fname)],
		key=os.path.getmtime)[-1]

if not os.path.isdir(htmldir) or \
        os.path.getmtime(newest_rst) > os.path.getmtime(newest_html):
    shoutout('sphinx html documentation out of date!')

# additional data files {{{1
docs_datafiles = []
for path, dirs, files in os.walk(htmldir):
    docs_datafiles += [(path, [os.path.join(path, fname) for fname in files])]

testsuites_datafiles = []
for path, dirs, files in os.walk(os.path.join(os.path.dirname(__file__), 'testsuites')):
    testsuites_datafiles += [(path, [os.path.join(path, fname) for fname in files])]

# setup(... {{{1
setup(
    name='kvarq',
    version=kvarq.VERSION,
    description='Targeted and Direct Variant Calling from FastQ Reads of Bacterial Genomes',
    long_description='KvarQ performs rapid in-silico genotyping for selected loci (e.g. phylogenetic SNPs, drug resistance mutations, repeats) in bacterial genome sequences in FastQ format. Mapping to a whole-genome reference sequence or de-novo assembly or the short reads is not necessary.',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Scientists, Clinicians',
        'License :: OSI Approved :: GPL License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        ],
    keywords='',
    author='Andreas Steiner',
    author_email='andreas.steiner@unibas.ch',
    url=kvarq.DOWNLOAD_URL,
    license='GPL v3',
    package_dir={'kvarq': 'kvarq'}, 
    packages = [
        'kvarq',
        'kvarq.gui'
        ],
    ext_modules=[engine_ext],
    test_suite = 'tests.automated',
    # (package_data does not work well with py2app/py2exe packaged applications)
    data_files = docs_datafiles + testsuites_datafiles + additional_datafiles,
    include_package_data=False, # from version control...
    app=[os.path.join('kvarq', 'gui', 'main.py')],
    entry_points={
        'console_scripts': [
             'kvarq = kvarq.cli:main'
        ],
    },
    **kwargs
#    install_requires=[],
#    tests_require=[],
)


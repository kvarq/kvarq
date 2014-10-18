"""
since KvarQ 0.12.2 testsuites are loaded in a two phase process

  1. testsuites are discovered using ``discover_testsuites``
  2. then testsuites are properly loaded using ``load_testsuites``
     or ``update_testsuites``

the :py:func:`kvarq.testsuites.load_testsuites` mechanism should
be used when a specified selection of testsuites is loaded, e.g.
for scanning. it simply selects testsuites from the pool of
discovered testsuites by their name.

the :py:func:`kvarq.testsuites.update_testsuites` mechanism should
be used when specific testsuites with a given version need to be
loaded, e.g. for exploring a ``.json`` file.
"""

import os.path, os, time
from os.path import expanduser

from distutils.version import StrictVersion

from util import get_root_path
from log import lo
import genes
from analyse import TestsuiteVersionConflictException


def add_testsuites_dir(testsuite_paths, base):

    if not os.path.isdir(base):
        return

    for subdir in os.listdir(base):

        if not os.path.isdir(os.path.join(base, subdir)) or (
                subdir[0] == '_' or subdir[0] == '.'):
            continue

        for fname in os.listdir(os.path.join(base, subdir)):

            if not fname.endswith('.py') or (
                    fname[0] == '_' or fname[0] == '.'):
                continue

            name = subdir + '/' + fname[:-3]
            path = os.path.join(base, subdir, fname)
            if name in testsuite_paths:
                lo.info('testsuite %s loaded from "%s"' % (name, path))
            else:
                lo.debug('testsuite %s loaded from "%s"' % (name, path))
            testsuite_paths[name] = path


def discover_testsuites(paths=[]):
    ''' returns dictionary mapping name to python file for all
        testsuites discovered in the usual places: kvarq root path,
        user home directory, current working directory,
        KVARQ_TESTSUITES environment variable, and any more paths
        specified paths as arguments -- later occurrences of the same
        testsuite override previous '''

    testsuite_paths = {}

    # 1) discover in root path
    root_base = os.path.abspath(os.path.join(get_root_path(), 'testsuites'))
    lo.debug('discovering testsuites in root path')
    add_testsuites_dir(testsuite_paths, root_base)

    # 2) discover from $HOME
    base = os.path.join(expanduser('~'), 'kvarq_testsuites')
    lo.debug('discovering testsuites in home directory')
    add_testsuites_dir(testsuite_paths, base)

    # 3) discover from CWD if not in root path
    cwd_base = os.path.abspath('testsuites')
    if cwd_base != root_base:
        lo.debug('discovering testsuites in current working directory')
        add_testsuites_dir(testsuite_paths, cwd_base)

    # 4) discover from KVARQ_TESTSUITES
    from_env = os.environ.get('KVARQ_TESTSUITES')
    if from_env:
        lo.debug('discovering testsuites in $KVARQ_TESTSUITES')
        for base in from_env.split(os.path.pathsep):
            add_testsuites_dir(testsuite_paths, base)

    # 5) explicitely specified paths
    for base in paths:
        if os.path.isdir(base):
            lo.debug('discovering testsuites in "%s"' % base)
            add_testsuites_dir(testsuite_paths, base)
        else:
            lo.warning('could not find directory "%s"' % base)

    return testsuite_paths


def load_testsuite(path):
    t0 = time.time()
    testsuite = genes.load_testsuite(path)
    lo.info('loaded testsuite from "%s" in %dms' % (
        path, int(1e3*(time.time() - t0))))
    return testsuite


def load_get_testsuite(testsuites, name, testsuite_paths):
    ''' name can be full name or only filename part '''

    if name in testsuites:
        return testsuites[name]

    if name in testsuite_paths:
        return load_testsuite(testsuite_paths[name])

    for fullname, path in testsuite_paths.items():
        if fullname.split('/')[-1] == name:

            lo.info('mapping testsuite "%s" to "%s"' % (name, fullname))

            if fullname in testsuites:
                return testsuites[fullname]

            return load_testsuite(testsuite_paths[fullname])

    return None


def load_testsuites(testsuite_paths, selection, raise_exception=False):
    ''' loads testsuites specified by their full name, group name,
        or directly specify testsuite by path of python file '''

    testsuites = {}

    groups = {}
    for name, path in testsuite_paths.items():
        parts = name.split('/')
        groups.setdefault(parts[0], {})[name] = path

    for name_or_path in selection:

        try:

            if (os.path.isfile(name_or_path)
                    and not name_or_path in testsuite_paths
                    and not name_or_path in groups):

                # load from path
                parts = name_or_path.split(os.path.sep)
                name = os.path.splitext(parts[-1])[0]
                if len(parts) > 1:
                    name = parts[-2] + '/' + name
                path = name_or_path

                testsuites[name] = load_testsuite(path)

            elif name_or_path in groups:

                # load group
                for name, path in groups[name_or_path].items():
                    testsuites[name] = load_testsuite(path)

            else:

                # load by full name
                name = name_or_path
                if name in testsuite_paths:
                    testsuites[name] = load_testsuite(testsuite_paths[name])
                else:
                    lo.warning('could not find testsuite "%s"' % name)

        except genes.TestsuiteLoadingException as e:
            if raise_exception:
                raise e
            lo.error('could not load testsuite from "%s" : %s' % (name, e))

    return testsuites


def update_testsuites(testsuites, names_versions, testsuite_paths):
    ''' updates the testsuites dictionary by adding specified testsuites
        if they are not yet loaded. raises TestsuiteLoadingException
        if a testsuite cannot be found and TestsuiteVersionConflictException

        :param testsuites: dictionary to update; same value can appear
            with different keys (see below)

        :param names_versions: array specifying testsuite name and
            version to be loaded; name can be full name (e.g.
            ``MTBC/resistance``) or partial name (e.g. ``resistance``,
            for backwards compatability) -- in any case, the testsuites
            dictionary will be updated with the name as specified in
            this argument

        :param testsuite_paths: dictionary as returned by
            :py:func:`kvarq.testsuites.discover_testsuites``

    '''

    namemap = dict([
            (name.split('/')[-1], name) for name in testsuite_paths
        ])

    for name, version in names_versions.items():

        testsuite = load_get_testsuite(testsuites, name, testsuite_paths)

        if testsuite is None:
            raise genes.TestsuiteLoadingException(
                    'could not find testsuite "%s"' % name)

        v = StrictVersion(version)
        tv = StrictVersion(testsuite.version)

        if tv < v or tv.version[0] != v.version[0]:
            raise TestsuiteVersionConflictException(
                    'incompatible versions testsuite "%s" : '
                    'expected %s found %s' % (name, v, tv))

        testsuites[name] = testsuite


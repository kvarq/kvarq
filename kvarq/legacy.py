
from distutils.version import StrictVersion
from collections import OrderedDict

from kvarq import VERSION
from kvarq import genes
from kvarq.log import lo


def convert_legacy_data(testsuites, data):
    '''
    :param testsuites: dictionary of :py:class:`kvarq.genes.Testsuite`
    :param data: dictionary as returned by :py:meth:`kvarq.analyse.Analyser.encode`,
        which can be of a previous version of KvarQ
    :returns: object as returned by :py:meth:`kvarq.analyse.Analyser.encode` of
        current version of KvarQ

    tries to convert data from older versions of KvarQ, raises
    :py:class:`kvarq.analyse.VersionConflictException` or
    :py:class:`kvarq.analyse.DataInconcistencyException` if data cannot be converted
    '''
    from kvarq.analyse import VersionConflictException, DataInconcistencyException

    kvarq_version = list(StrictVersion(VERSION).version)
    version = list(StrictVersion(data['info']['version']).version)

    if version[1] < 10:
        raise VersionConflictException('cannot load files v<0.10')

    # convert tests -> coverages
    if version[0] == 0 and version[1] == 10:

        # load data
        templates_by_testname = dict(reduce(lambda x, y: x + y, [[
                (str(test), test.template) for test in testsuite.tests
            ] for testsuite in testsuites.values()]))

        coverages_by_testname = dict(reduce(lambda x, y: x + y,
                [data_testsuite.items() for data_testsuite in data['testsuites'].values()]
            ))

        # convert test-nr -> coverage-nr
        nrmap = []
        coverages = OrderedDict()

        for i, testname in enumerate(data['tests']):

            if not testname in templates_by_testname:
                lo.info('json contains additional test "%s"; discarding.' % testname)
                continue

            templatename = str(templates_by_testname[testname])
            coverage = coverages_by_testname[testname]

            if templatename in coverages:
                assert coverages[templatename] == coverage, DataInconcistencyException(
                        'found contradicting coverages for template "%s" : "%s" / "%s"' %
                        templatename, (coverages[templatename], coverage))
            else:
                coverages[templatename] = coverage
                nrmap.append(i)

        # save to data
        data['coverages'] = [(k, v) for k, v in coverages.items()]
        lo.debug('mapping "nseqhits", "nseqbasehits" : (%d) %s' % (len(nrmap), str(nrmap)))
        for key in ['nseqhits', 'nseqbasehits']:
            data['stats'][key] = [
                    data['stats'][key][nrmap[coveragenr]] 
                    for coveragenr in range(len(coverages)) # forward
                ] + [
                    data['stats'][key][nrmap[coveragenr] + len(data['tests'])] 
                    for coveragenr in range(len(coverages)) # reverse
                ]

        # clean up
        del data['testsuites']
        del data['tests']
        version[1] += 1

    assert version[0] == kvarq_version[0] and version[1] == kvarq_version[1], \
            VersionConflictException('could not elevate version more than to "%d.%d"' %
                    (version[0], version[1]))

    return data

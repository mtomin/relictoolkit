from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import next
from future import standard_library
import MDAnalysis
import relictoolkit.calculation as c
import os
from configparser import ConfigParser
standard_library.install_aliases()


def test_add_lj_parameters():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    c.add_lj_parameters(system)
    assert system.atoms[0].lj_energy == 3.66


def test_perform_analysis():
    test_config = ConfigParser()
    test_config['files'] = {
        'topology': os.path.dirname(__file__) + '/data/testtop.prmtop',
        'trajectories': os.path.dirname(__file__) + '/data/testtraj.xcrd',
        'filetype': 'mdcrd'
    }
    test_config['parameters'] = {
        'mask1': 'resid 3',
        'mask2': 'resid 4',
        'stride': 1,
        'ncores': 1,
        'dt': 2
    }

    with open('testrun_config.ini', 'w+') as f:
        test_config.write(f)

    c.perform_analysis('testrun_config.ini')
    logfile = open('relic_logfile.log', 'r+')
    log_line = logfile.readline()
    assert log_line.split()[0] == 'Calculation'
    next(logfile)
    next(logfile)
    log_line = logfile.readline()
    assert log_line.split()[0] == 'Topology:'
    logfile.close()

    with open('output.dat', 'r+') as output:
        next(output)
        next(output)
        output_line = output.readline()
        assert output_line.split()[2] == '-3.22555'

    os.remove('output.dat')
    os.remove('testrun_config.ini')

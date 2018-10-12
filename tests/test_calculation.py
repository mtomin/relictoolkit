from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
import MDAnalysis
import relictoolkit.calculation as c
import os
standard_library.install_aliases()


def test_add_lj_parameters():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    c.add_lj_parameters(system)
    assert system.atoms[0].lj_energy == 3.66


def test_process_trajectory():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    mask1 = 'resid 3'
    mask2 = 'resid 4'
    try:
        os.symlink(os.path.dirname(__file__) + '/data/testtraj.xcrd', os.path.dirname(__file__) + '/data/testtraj.mdcrd')
    except FileExistsError:
        os.remove(os.path.dirname(__file__) + '/data/testtraj.mdcrd')
        os.symlink(os.path.dirname(__file__) + '/data/testtraj.xcrd', os.path.dirname(__file__) + '/data/testtraj.mdcrd')

    if os.path.isfile('relic_logfile.log'):
        os.remove('relic_logfile.log')

    c.process_trajectory(topology, [os.path.dirname(__file__) + '/data/testtraj.mdcrd'], 2, 1, 1, mask1, mask2, 0)
    logfile = open('relic_logfile.log', 'r+')
    log_line = logfile.readline()
    print(log_line)
    assert log_line == 'Core 0 assigned frames 0 to 1\n'
    log_line = logfile.readline()
    print(log_line)
    assert log_line == 'Processing trajectory segment 0 frame 1 of 1\n'
    logfile.close()
    os.remove('relic_logfile.log')
    with open('output.dat_0', 'r+') as output:
        next(output)
        next(output)
        output_line = output.readline()
        assert output_line.split()[2] == '-3.22555'
    os.remove('output.dat_0')

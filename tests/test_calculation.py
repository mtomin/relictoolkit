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
import mock
from configparser import ConfigParser
standard_library.install_aliases()


def test_add_lj_parameters():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    c.add_lj_parameters(system)
    assert system.atoms[0].lj_energy == 3.66
    assert system.atoms[0].vdw_radius == 0.069


def test_timeit():

    @c.timeit
    def test_function():
        pass

    if os.path.isfile('relic_logfile.log'):
        os.remove('relic_logfile.log')

    test_function()
    with open('relic_logfile.log', 'r+') as logfile:
        next(logfile)
        line = logfile.readline()
    assert line == 'Calculation finished successfully\n'
    os.remove('relic_logfile.log')


def do_nothing(*args):
    pass


@mock.patch('relictoolkit.calculation.add_lj_parameters', side_effect=do_nothing)
@mock.patch('relictoolkit.utils.process_frame', side_effect=do_nothing)
@mock.patch('relictoolkit.utils.load_partial_traj')
@mock.patch('relictoolkit.utils.load_from_config')
def test_process_trajectory(mock_load_config, mock_load_partial_traj, mock_process_frame, mock_add_lj_parameters):
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

    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    mock_load_config.return_value = ['output.dat']
    mock_load_partial_traj.return_value = {'traj': system.trajectory,
                                           'startframe': 0,
                                           'endframe': 1}
    c.process_trajectory(topology, [os.path.dirname(__file__) + '/data/testtraj.mdcrd'], 2, 1, 1, mask1, mask2, 0)
    logfile = open('relic_logfile.log', 'r+')
    log_line = logfile.readline()

    assert log_line == 'Processing trajectory segment 0 frame 1 of 1\n'
    logfile.close()
    os.remove('relic_logfile.log')
    with open('output.dat_0', 'r+') as output:
        output_line = output.readline()

    assert output_line == 'Frame #    Residue #     Eelec      Evdw       Etotal     Timestep: 2\n'
    os.remove('output.dat_0')
    os.remove(os.path.dirname(__file__) + '/data/testtraj.mdcrd')


@mock.patch('relictoolkit.utils.write_logfile_header', side_effect=do_nothing())
@mock.patch('multiprocessing.Pool', side_effect=do_nothing())
def test_perform_analysis(*args):
    # Create a test config file
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
    os.remove('testrun_config.ini')

    assert os.path.isfile('output.dat')
    assert not os.path.isfile('output.dat_0')
    with open('output.dat') as outfile:
        next(outfile)
        next(outfile)
        line = outfile.readline()
        assert line.split()[2] == '-1.26703'
    os.remove('output.dat')

test_process_trajectory()
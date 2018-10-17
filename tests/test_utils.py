from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from future import standard_library
import relictoolkit.utils as u
import relictoolkit.calculation as p
import MDAnalysis
from configparser import ConfigParser
import os
import mock
import shutil
standard_library.install_aliases()


def test_vdw_energy():
    topology = os.path.dirname(__file__)+'/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__)+'/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    p.add_lj_parameters(system)
    energy = u.vdw_energy(system.atoms[0], system.atoms[1], system.dimensions)
    assert energy == -7.1872471955340208e-07


def test_electrostatic_energy():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    energy = u.electrostatic_energy(system.atoms[0], system.atoms[1], system.dimensions)
    assert energy == 0.82850458551545669


def test_load_from_config():
    test_config = os.path.dirname(__file__)+'/data/test_config.ini'
    assert u.load_from_config('files', 'filetype', test_config) == ['mdcrd', False]
    assert u.load_from_config('files', 'output', test_config) == ['/testfolder/testoutput.out', False]
    assert u.load_from_config('files', 'topology', test_config) == ['/testfolder/testtopology.prmtop', False]
    assert u.load_from_config('files', 'trajectories', test_config) == ['/testfolder/testtrajectory.xcrd', False]
    assert u.load_from_config('parameters', 'indi', test_config) == ['1', False]
    assert u.load_from_config('parameters', 'cutoff', test_config) == ['1', False]
    assert u.load_from_config('parameters', 'vdw', test_config) == ['True', False]
    assert u.load_from_config('parameters', 'mask1', test_config) == ['resid 1', False]
    assert u.load_from_config('parameters', 'ncores', test_config) == ['1', False]
    assert u.load_from_config('parameters', 'stride', test_config) == ['1', False]
    assert u.load_from_config('parameters', 'mask2', test_config) == ['resid 2', False]


def test_load_config_defaults():
    defaults_config = os.path.dirname(__file__)+'/data/test_config_defaults.ini'
    assert u.load_from_config('files', 'output', defaults_config) == ['output.dat', True]
    assert u.load_from_config('parameters', 'indi', defaults_config) == ['4', True]
    assert u.load_from_config('parameters', 'cutoff', defaults_config) == ['5', True]
    assert u.load_from_config('parameters', 'vdw', defaults_config) == ['False', True]
    assert u.load_from_config('parameters', 'ncores', defaults_config) == ['1', True]
    assert u.load_from_config('parameters', 'stride', defaults_config) == ['1000', True]


def test_load_plot_config():
    test_config_filename = os.path.dirname(__file__)+'/data/test_config_plot.ini'
    config = ConfigParser()
    config.read(test_config_filename)
    config.set('files', 'datafile', os.path.dirname(__file__)+'/data/test_output.out')
    with open(test_config_filename, 'w+') as f:
        config.write(f)

    assert u.load_from_plot_config('parameters', 'plot_type', test_config_filename) == 'averages'
    assert u.load_from_plot_config('parameters', 'startframe', test_config_filename) == '2'
    assert u.load_from_plot_config('parameters', 'starting_residue', test_config_filename) == '3'
    assert u.load_from_plot_config('parameters', 'endframe', test_config_filename) == '4'
    assert u.load_from_plot_config('parameters', 'end_residue', test_config_filename) == '5'

    config.set('files', 'datafile', '/data/test_output.out')
    with open(test_config_filename, 'w+') as f:
        config.write(f)


def test_load_plot_config_defaults():
    test_config_filename = os.path.dirname(__file__)+'/data/test_config_plot_defaults.ini'
    config = ConfigParser()
    config.read(test_config_filename)
    config.set('files', 'datafile', os.path.dirname(__file__)+'/data/test_output.out')
    with open(test_config_filename, 'w+') as f:
        config.write(f)

    assert u.load_from_plot_config('parameters', 'plot_type', test_config_filename) == 'time'
    assert u.load_from_plot_config('parameters', 'startframe', test_config_filename) == 0
    assert u.load_from_plot_config('parameters', 'starting_residue', test_config_filename) == 1
    assert u.load_from_plot_config('parameters', 'endframe', test_config_filename) == 4500
    assert u.load_from_plot_config('parameters', 'end_residue', test_config_filename) == 692

    config.set('files', 'datafile', '/data/test_output.out')
    with open(test_config_filename, 'w+') as f:
        config.write(f)


def test_check_params():
    shutil.copyfile(os.path.dirname(__file__) + '/data/test_config.ini',
                    os.path.dirname(__file__) + '/data/test_config.ini_bak')
    try:
        config = ConfigParser()
        config_filename = os.path.dirname(__file__) + '/data/test_config.ini'
        config.read(config_filename)
        assert u.check_params(config) == 'Topology file missing!'
        config.set('files', 'topology', os.path.dirname(__file__) + '/data/testtop.prmtop')
        config.set('files', 'trajectories', os.path.dirname(__file__) + '/data/testtraj.xcrd')
        assert u.check_params(config, config_filename) == ''
        config.set('parameters', 'mask1', 'klj')
        assert u.check_params(config, config_filename) == 'Mask1 error!'
        config.set('parameters', 'mask1', 'resid 1 and resid 2')
        assert u.check_params(config, config_filename) == 'Mask1 contains no atoms!'
        config.set('parameters', 'mask1', 'resid 1')
        config.set('parameters', 'mask2', 'klj')
        assert u.check_params(config, config_filename) == 'Mask2 error!'
        config.set('parameters', 'mask2', 'resid 1 and resid 2')
        assert u.check_params(config, config_filename) == 'Mask2 contains no atoms!'
        config.set('parameters', 'mask2', 'resid 2')
        config.set('parameters', 'ncores', '300')
        with open(config_filename, 'w+') as f:
            config.write(f)
        assert u.check_params(config, config_filename).split('(')[0] == \
               'Number of cores specified higher than available number of cores '

        assert u.check_params(config, config_filename).split('(')[0] == 'Cutoff must be a number!'
        config.set('parameters', 'ncores', '1')
        config.set('parameters', 'cutoff', 'asd')
        with open(config_filename, 'w+') as f:
            config.write(f)

        shutil.move(os.path.dirname(__file__) + '/data/test_config.ini_bak',
                    os.path.dirname(__file__) + '/data/test_config.ini')
    except:
        shutil.move(os.path.dirname(__file__) + '/data/test_config.ini_bak',
                    os.path.dirname(__file__) + '/data/test_config.ini')


def test_check_plot_params():
    config = ConfigParser()
    config_filename = os.path.dirname(__file__) + '/data/test_config_plot.ini'
    config.read(config_filename)
    config.set('files', 'datafile', os.path.dirname(__file__) + '/data/test_output.out')
    with open(config_filename, 'w+') as configfile:
        config.write(configfile)
    assert u.check_plot_params(config, config_filename) == ''
    config.set('parameters', 'startframe', 'asdf')
    with open(config_filename, 'w+') as f:
        config.write(f)
    assert u.check_plot_params(config, config_filename) == 'Starting frame must be an integer!'
    config.set('parameters', 'startframe', '2')
    with open(config_filename, 'w+') as f:
        config.write(f)
    config.set('files', 'datafile', 'teststr.ini')
    assert u.check_plot_params(config, config_filename) == 'File not found!'

    config.set('files', 'datafile', '/data/test_output.out')
    with open(config_filename, 'w+') as configfile:
        config.write(configfile)


@mock.patch('relictoolkit.utils.electrostatic_energy')
def test_interdomain_interactions(electrostatic_energy_mock):
    electrostatic_energy_mock.return_value = -1
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    domain1 = system.select_atoms('resid 1')
    domain2 = system.select_atoms('resid 3')
    interaction = u.interdomain_interactions(domain1, domain2, 1)
    assert interaction[0].split()[2] == '-285.00000'


def test_read_uff_parameters():
    uff_parameters = u.read_uff_parameters()
    assert uff_parameters['y'] == [3.345, 0.072]


def test_write_logfile_header():
    config = ConfigParser()
    config.read(os.path.dirname(__file__) + '/data/test_config.ini')
    with open('testlogfile', 'w+') as testlogfile:
        u.write_logfile_header(config, testlogfile)

    with open('testlogfile', 'r+') as testlogfile:
        assert testlogfile.readline()[0:19] == 'Calculation started'


def test_load_partial_traj():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    result = u.load_partial_traj(system, 1, 1, 0)
    assert result['traj'][0].dimensions[3] == 90.0
    assert result['startframe'] == 0


@mock.patch('relictoolkit.utils.interdomain_interactions')
def test_process_frame(interdomain_interactions_mock):
    interdomain_interactions_mock.return_value = ['0 1 2 0 2']
    with open('test_frameprocess.dat', 'w+') as testoutput:
        u.process_frame('resid 1', 'resid2', testoutput, 1)

    with open('test_frameprocess.dat', 'r+') as testoutput:
        line = testoutput.readline()
    assert line == '0 1 2 0 2\n'
    os.remove('test_frameprocess.dat')


def test_tail():
    with open('testtail.dat', 'w+') as testtail:
        for i in range(0, 5):
            print('%s\n' % i, file=testtail)
        print('the end', file=testtail)
        assert u.tail(testtail)[0] == 'the end\n'
        os.remove('testtail.dat')

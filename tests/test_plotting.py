from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import next
from builtins import int
from future import standard_library
from configparser import ConfigParser
import os
import relictoolkit.plotting as p
import mock
import shutil
standard_library.install_aliases()


def test_read_datafile():
    with open(os.path.dirname(__file__) + '/data/test_output.out') as data:
        next(data)
        dt = int(data.readline().split()[1])
        datapoints = p.read_datafile(data, 0, 4500, 1, 692, 'time', dt)
    assert datapoints == {'z': [0.0, None, -38.58544], 'y': [1.0, None, 692.0], 'step': 4500, 'x': [0.0, None, 2.25]}

    with open(os.path.dirname(__file__) + '/data/test_output.out') as data:
        next(data)
        dt = int(data.readline().split()[1])
        datapoints = p.read_datafile(data, 0, 4500, 1, 692, 'averages', dt)
    assert datapoints == {'x': [None], 'z': [0.0, None, -38.58544], 'y': [1.0, None, 692.0], 'step': 4500}


def test_generate_layout():
    test_layout = p.generate_layout('time')
    assert test_layout['title'] == 'Residue interactions'
    assert test_layout['scene']['zaxis']['title'] == 'E/kJmol<sup>-1</sup>'
    assert test_layout['scene']['xaxis']['title'] == 'Time/ns'

    test_layout = p.generate_layout('frame_number')
    assert test_layout['scene']['xaxis']['title'] == 'Frame'

    test_layout = p.generate_layout('averages')
    assert test_layout['scene']['xaxis']['title'] == 'Residue'


def test_calculate_averages():
    datapoints = {'z': [0.0, -38.58544], 'y': [1.0, 692.0], 'step': 4500, 'x': [0.0, 2.25]}
    test_averages = p.calculate_averages(os.path.dirname(__file__) + '/data/test_output.out', datapoints, 1, 1, 1, 692)
    assert len(test_averages['residues']) == 692
    assert test_averages['residue_energies'][1] == -38.58544


@mock.patch('relictoolkit.plotting.read_datafile')
@mock.patch('relictoolkit.plotting.calculate_averages')
def test_generate_figure_data_plotly(calculate_averages_mock, read_datafile_mock):

    shutil.copy(os.path.dirname(__file__) + '/data/test_config_plot.ini', os.path.dirname(__file__) +
                '/data/test_config_plot.ini_bak')
    try:
        # Set path in confg_plot.ini
        test_config_filename = os.path.dirname(__file__) + '/data/test_config_plot.ini'
        config = ConfigParser()
        config.read(test_config_filename)
        config.set('files', 'datafile', os.path.dirname(__file__) + '/data/test_output.out')
        with open(test_config_filename, 'w+') as f:
            config.write(f)

        read_datafile_mock.return_value = {
            'z': [0.0, None, -38.58544], 'y': [1.0, None, 692.0], 'step': 4500, 'x': [0.0, None, 2.25]
        }
        calculate_averages_mock.return_value = {
            'residues': [0.0, 1.0, 2.0],
            'residue_energies': [0.0, 5.0, 10.0]
        }
        testfigure = p.generate_figure_data_plotly(os.path.dirname(__file__) + '/data/test_config_plot.ini')

        assert testfigure['layout']['title'] == 'Average residue energies'
        assert testfigure['data'][0]['x'] == (0.0, 1.0, 2.0)
        assert testfigure['data'][0]['y'] == (0.0, 5.0, 10.0)

        config.set('parameters', 'plot_type', 'time')
        with open(test_config_filename, 'w+') as f:
            config.write(f)

        testfigure = p.generate_figure_data_plotly(os.path.dirname(__file__) + '/data/test_config_plot.ini')
        assert testfigure['data']['z'] == (0.0, None, -38.58544)

        shutil.copy(os.path.dirname(__file__) + '/data/test_config_plot.ini_bak', os.path.dirname(__file__) +
                    '/data/test_config_plot.ini')
    except:
        shutil.copy(os.path.dirname(__file__) + '/data/test_config_plot.ini_bak', os.path.dirname(__file__) +
                    '/data/test_config_plot.ini')


@mock.patch('relictoolkit.plotting.read_datafile')
@mock.patch('relictoolkit.plotting.calculate_averages')
def test_generate_figure_data_mplt(calculate_averages_mock, read_datafile_mock):

    test_config_filename = os.path.dirname(__file__) + '/data/test_config_plot.ini'
    config = ConfigParser()
    config.read(test_config_filename)
    config.set('files', 'datafile', os.path.dirname(__file__) + '/data/test_output.out')
    with open(test_config_filename, 'w+') as f:
        config.write(f)

    calculate_averages_mock.return_value = {
        'residues': [0.0, 1.0, 2.0],
        'residue_energies': [0.0, 5.0, 10.0]
    }

    testfigure = p.generate_figure_data_mplt(test_config_filename)
    testfigure.show()
    assert testfigure.gca().lines[0].get_xdata()[1] == 1.0
    assert testfigure.gca().lines[0].get_ydata()[2] == 10.0

test_generate_figure_data_plotly()
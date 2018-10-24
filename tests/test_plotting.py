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
        dt = int(data.readline().split()[-1])
        datapoints = p.read_datafile(data, 0, 4500, 1, 692, 'time', dt)

    assert datapoints == {'z': [0.0, None, -38.58544], 'y': [1.0, None, 692.0], 'step': 4500, 'x': [0.0, None, 9.0]}

    with open(os.path.dirname(__file__) + '/data/test_output.out') as data:
        next(data)
        dt = int(data.readline().split()[1])
        datapoints = p.read_datafile(data, 0, 4500, 1, 692, 'averages', dt)

    assert datapoints == {'x': [None], 'z': [None, -38.58544], 'y': [None, 692.0], 'step': 4500}

    with open(os.path.dirname(__file__) + '/data/test_output.out') as data:
        dt = int(data.readline().split()[-1])
        datapoints = p.read_datafile(data, 0, 4500, 1, 692, 'frame_number', dt)

    assert datapoints =={'z': [0.0, None, -38.58544], 'x': [0.0, None, 4500.0], 'step': 4500, 'y': [1.0, None, 692.0]}


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
    assert len(test_averages['residues']) == 691
    assert test_averages['residue_energies'][1] == -38.58544


@mock.patch('relictoolkit.plotting.read_datafile')
@mock.patch('relictoolkit.plotting.calculate_averages')
def test_generate_figure_data_plotly(calculate_averages_mock, read_datafile_mock):

    shutil.copy(os.path.dirname(__file__) + '/data/test_config_plot.ini', os.path.dirname(__file__) +
                '/data/test_config_plot_temp.ini')

    # Set path in confg_plot.ini
    test_config_filename = os.path.dirname(__file__) + '/data/test_config_plot_temp.ini'
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
    testfigure = p.generate_figure_data_plotly(os.path.dirname(__file__) + '/data/test_config_plot_temp.ini')

    assert testfigure['layout']['title'] == 'Average residue energies'
    assert testfigure['data'][0]['x'] == (0.0, 1.0, 2.0)
    assert testfigure['data'][0]['y'] == (0.0, 5.0, 10.0)

    config.set('parameters', 'plot_type', 'time')
    with open(test_config_filename, 'w+') as f:
        config.write(f)

    testfigure = p.generate_figure_data_plotly(os.path.dirname(__file__) + '/data/test_config_plot_temp.ini')
    assert testfigure['data'][0]['z'] == (0.0, None, -38.58544)

    config.set('parameters', 'plot_type', 'frame_number')
    with open(test_config_filename, 'w+') as f:
        config.write(f)
    assert testfigure['data'][0]['z'] == (0.0, None, -38.58544)
    assert testfigure['data'][0]['x'] == (0.0, None, 2.25)

    os.remove(os.path.dirname(__file__) + '/data/test_config_plot_temp.ini')


@mock.patch('relictoolkit.plotting.read_datafile')
@mock.patch('relictoolkit.plotting.calculate_averages')
def test_generate_figure_data_mplt(calculate_averages_mock, read_datafile_mock):

    shutil.copy(os.path.dirname(__file__) + '/data/test_config_plot.ini', os.path.dirname(__file__) +
                '/data/test_config_plot_temp.ini')

    test_config_filename = os.path.dirname(__file__) + '/data/test_config_plot_temp.ini'
    config = ConfigParser()
    config.read(test_config_filename)
    config.set('files', 'datafile', os.path.dirname(__file__) + '/data/test_output.out')
    with open(test_config_filename, 'w+') as f:
        config.write(f)

    read_datafile_mock.return_value = {
        'z': [0.0, 15.0, None, 1.0, 16.0], 'y': [0.0, 1.0, None, 0.0, 1.0], 'step': 4500, 'x': [0.0, 0.0, None, 0.1, 0.1]
    }

    calculate_averages_mock.return_value = {
        'residues': [0.0, 1.0, 2.0],
        'residue_energies': [0.0, 5.0, 10.0]
    }

    testfigure = p.generate_figure_data_mplt(test_config_filename)
    assert testfigure.gca().lines[0].get_xdata()[1] == 1.0
    assert testfigure.gca().lines[0].get_ydata()[2] == 10.0

    config.set('parameters', 'plot_type', 'time')
    with open(test_config_filename, 'w+') as f:
        config.write(f)

    testfigure = p.generate_figure_data_mplt(test_config_filename)
    assert testfigure.get_axes()[0].get_children()[9].get_text() == 'Residue interactions'

    os.remove(os.path.dirname(__file__) + '/data/test_config_plot_temp.ini')


def do_nothing(*args, **kwargs):
    pass


@mock.patch('relictoolkit.utils.load_from_plot_config', side_effect=['averages', 'False', 'time', 'True'])
@mock.patch('plotly.offline.plot', side_effect=do_nothing)
@mock.patch('relictoolkit.plotting.generate_figure_data_mplt')
@mock.patch('matplotlib.pyplot.figure', 'show')
@mock.patch('relictoolkit.plotting.generate_figure_data_plotly', side_effect=do_nothing)
def test_plotting_main(generate_data_plotly, pyplot_figure_show_mock, generate_data_mplt_mock, load_from_plot_config_mock):
    p.main()
    assert pyplot_figure_show_mock.call_args[0][0] == 'config_plot.ini'
    assert generate_data_mplt_mock.call_args is None
    assert load_from_plot_config_mock.call_args[0] == ('parameters', 'interactive', 'config_plot.ini')

    p.main()
    assert pyplot_figure_show_mock.call_args[0][0] == 'config_plot.ini'
    assert load_from_plot_config_mock.call_args[0] == ('parameters', 'interactive', 'config_plot.ini')

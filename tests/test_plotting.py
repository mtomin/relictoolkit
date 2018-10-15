from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import next
from builtins import int
from future import standard_library
import os
import relictoolkit.plotting as p
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


test_read_datafile()
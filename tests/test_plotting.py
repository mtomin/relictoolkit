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

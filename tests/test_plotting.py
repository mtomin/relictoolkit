import os
import relictoolkit.plotting as p


def test_read_datafile():
    with open(os.path.dirname(__file__) + '/data/test_output.out') as data:
        next(data)
        dt = int(data.readline().split()[1])
        datapoints = p.read_datafile(data, 0, 4500, 1, 692, 'time', dt)
    assert datapoints == {'z': [0.0, None, -38.58544], 'y': [1.0, None, 692.0], 'step': 4500, 'x': [0.0, None, 2.25]}

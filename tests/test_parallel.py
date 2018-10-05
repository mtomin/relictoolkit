from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
import MDAnalysis
import relictoolkit.calculation as p
import os
standard_library.install_aliases()


def test_add_lj_parameters():
    topology = os.path.dirname(__file__) + '/data/testtop.prmtop'
    trajectory = os.path.dirname(__file__) + '/data/testtraj.xcrd'
    system = MDAnalysis.Universe(topology, trajectory, format='mdcrd')
    p.add_lj_parameters(system)
    assert system.atoms[0].lj_energy == 3.66

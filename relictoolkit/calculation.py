from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from builtins import range
from builtins import int
from builtins import open
from builtins import str
from future import standard_library
from configparser import ConfigParser
import time
import numpy as np
import relictoolkit.utils as u
import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import AtomAttr
from multiprocessing import Pool
from functools import partial
import glob
import os
import shutil
import re
standard_library.install_aliases()


class LJenergy(AtomAttr):
    """Define Lennard-Jones parameters as a Universe class attribute
    """

    attrname = 'lj_energies'
    singular = 'lj_energy'


class Vdwradius(AtomAttr):
    """Define van der Waals radii as a Universe class attribute
    """

    attrname = 'vdw_radii'
    singular = 'vdw_radius'


def add_lj_parameters(system):
    """Append vdW radii and LJ parameters to all the atoms

    Parameters
    ----------
    system: MDAnalysis.Universe
        Loaded system
    """

    uff_parameters = u.read_uff_parameters()
    energy_parameters = list()
    vdw_radii = list()

    # TODO: Test more formats
    for atom in system.atoms:
        energy_parameters.append(uff_parameters[atom.element.lower()][0])
        vdw_radii.append(uff_parameters[atom.element.lower()][1])

    system.add_TopologyAttr('lj_energy', values=np.array(energy_parameters))
    system.add_TopologyAttr('vdw_radius', values=np.array(vdw_radii))


def process_trajectory(topology, trajlist, dt, step, ncores, mask1, mask2, core):
    """Iterates through partial trajectory frames and calculates interactions. Supports parallel processing by loading
    ncores trajectory segments.

    Parameters
    ----------
        topology : str
            Topology file (or standalone format)
        trajlist : list of str
            Filenames of trajectory files
        dt: int
            Timestep used in MD simulations
        step : int
            Step between two sampled trajectory frames
        ncores : int
            Number of CPU cores
        mask1 : str
            Selection containing the atoms in the first segment of interest
        mask2 : str
            Selection containing the atoms in the second segment of interest
        core : int
            CPU core used
    """

    input_system = mda.Universe(topology, trajlist, dt=dt)
    partial_traj = u.load_partial_traj(input_system, step, ncores, core)
    domain1 = input_system.select_atoms(mask1, updating=True)
    domain2 = input_system.select_atoms(mask2, updating=True)
    calculate_vdw = u.load_from_config('parameters', 'vdw') == 'True'
    if calculate_vdw:
        add_lj_parameters(input_system)
    else:
        input_system.add_TopologyAttr('lj_energy', values=[0]*len(input_system.atoms))
        input_system.add_TopologyAttr('vdw_radius', values=[0] * len(input_system.atoms))

    output = open(u.load_from_config('files', 'output')[0] + '_' + str(core), 'w+')

    # Print output header
    if core == 0:
        print('{:<9s}  {:<12s}  {:<10s} {:<10s} {:<10s} {:<10s}'.format(
            'Frame #', 'Residue #', 'Eelec', 'Evdw', 'Etotal', 'Timestep: %s' % dt), file=output)

    for structure in partial_traj['traj']:

        # Write current progress in logfile
        with open('relic_logfile.log', 'a+') as logfile:
            logfile.write(
                'Processing trajectory segment %s frame %s of %s\n' % (core, structure.frame//step -
                partial_traj['startframe']//step + 1, (partial_traj['endframe'] -
                                                       partial_traj['startframe'] - 1) // step + 1))

        # Calculate interactions for frame
        u.process_frame(domain1, domain2, output, structure.frame)
    output.close()


def perform_analysis(config_filename='config.ini'):
    """Loads settings from the config file, performs trajectory pre-processing and processing.
    """

    # Load input from config file
    config = ConfigParser()
    config.read(config_filename)
    trajectories = config.get('files', 'trajectories').split(', ')
    topology = config.get('files', 'topology')

    extension = u.load_from_config('files', 'filetype', config_filename)[0]
    stride = int(u.load_from_config('parameters', 'stride', config_filename)[0])
    mask1 = config.get('parameters', 'mask1')
    mask2 = config.get('parameters', 'mask2')
    ncores = int(u.load_from_config('parameters', 'ncores', config_filename)[0])
    dt = int(u.load_from_config('parameters', 'dt', config_filename)[0])

    # Check if the input is a standalone format
    standalone_format = True

    try:
        whole_traj = mda.Universe(topology, format=extension)
        testvar = len(whole_traj.trajectory)
    except (TypeError, ValueError) as error:
        standalone_format = False

    if standalone_format:
        trajectories = [topology]

    with open('relic_logfile.log', 'w+') as logfile:
        u.write_logfile_header(config, logfile)

    # Create temp symlinks for format issues
    for index, traj in enumerate(trajectories):
        temp_symlink_filename = 'temp_symlink_' + os.path.split(traj)[1] + '.' + extension

        try:
            os.symlink(traj, temp_symlink_filename)
        except FileExistsError:
            os.remove(temp_symlink_filename)
            os.symlink(traj, temp_symlink_filename)
        trajectories[index] = temp_symlink_filename

    with open('relic_logfile.log', 'a+') as logfile:
        logfile.write('Calculating interactions...\n')

    # Parallelization
    process_pool = Pool(processes=ncores)
    func = partial(process_trajectory, topology, trajectories, dt, stride, ncores, mask1, mask2)
    process_pool.map(func, range(ncores))
    process_pool.close()
    process_pool.join()

    for symlink in glob.glob('temp_symlink_*'):
        os.remove(symlink)

    # Merge partial outputs
    with open('relic_logfile.log', 'a+') as logfile:
        logfile.write('Merging outputs...\n')
    output = u.load_from_config('files', 'output', config_filename)[0]
    partial_outputs = glob.glob(output + '_*')

    # Sort files to avoid getting _1,12,13,..,19,2,21,..
    partial_outputs.sort(key=lambda x: int(re.search(r'(_\d*$)', x).group()[1:len(re.search(r'(_\d*$)', x).group())]))
    with open(output, 'wb') as final_output:
        for partial_output in partial_outputs:
            with open(partial_output, 'rb') as current_output:
                shutil.copyfileobj(current_output, final_output)
            os.remove(partial_output)


def timeit(program):
    """Time the function and append the runtime to logfile
    """
    def wrapper(*args):
        start_time = time.time()
        program(*args)
        logfile = open('relic_logfile.log', 'a+')
        logfile.write('\nCalculation finished successfully\n')
        logfile.write("Total calculation time: %s seconds" % (time.time() - start_time))
        logfile.close()
    return wrapper


@timeit
def main(config_filename='config.ini'):
    """Time the calculation and finalize the log file.
    """
    perform_analysis(config_filename)


if __name__ == "__main__":
    main()

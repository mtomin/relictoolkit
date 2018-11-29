from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import open
from builtins import int
from builtins import range
from future import standard_library
from configparser import ConfigParser
import numpy as np
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from math import pi
import datetime
import MDAnalysis
import multiprocessing
import os
import pkg_resources
standard_library.install_aliases()


def vdw_energy(atom1, atom2, box_dimensions):
    """Calculate van der Waals energy.
    Parameters
    ----------
    atom1: MDAnalysis.core.groups.Atom
        First atom
    atom2: MDAnalysis.core.groups.Atom
        Second atom
    box_dimensions: numpy.ndarray
        Dimensions of the system box (for PBC imaging)

    Returns
    -------
    eelec : float
        Van der Waals energy in kJ/mol
    """

    mixed_lj_term = np.sqrt(atom1.lj_energy*atom2.lj_energy)
    vdw_eq_radius = np.sqrt(atom1.vdw_radius*atom2.vdw_radius)
    r = MDAnalysis.lib.distances.calc_bonds(atom1.position.reshape(-1, 3),
                                            atom2.position.reshape(-1, 3), box_dimensions)[0]
    evdw = 4.184 * mixed_lj_term * (-2*(np.power(vdw_eq_radius/r, 6) + np.power(vdw_eq_radius/r, 12)))
    return evdw


def check_params(settings, filename='config.ini'):
    """Check validity of the parameters in the config.ini file.
    Parameters
    ----------
    settings: ConfigParser
        Object containing data located in config.ini
    filename: str
        Name of the config file
    """
    error_message = ''
    topology = settings.get('files', 'topology')
    trajectories = settings.get('files', 'trajectories').split(', ')
    extension = settings.get('files', 'filetype')
    mask1 = settings.get('parameters', 'mask1')
    mask2 = settings.get('parameters', 'mask2')
    cutoff = settings.get('parameters', 'cutoff')
    calculate_vdw = load_from_config('parameters', 'vdw')[0] == 'True'

    try:
        stride = int(load_from_config('parameters', 'stride', filename)[0])
    except ValueError:
        error_message = 'Stride must be an integer!'

    try:
        dt = int(load_from_config('parameters', 'dt', filename)[0])
    except ValueError:
        error_message = 'Timestep must be an integer!'

    try:
        ncores = int(load_from_config('parameters', 'ncores', filename)[0])
    except ValueError:
        error_message = 'Number of cores must be an integer!'

    try:
        whole_traj = MDAnalysis.Universe(topology)
    except (KeyError, ValueError) as error:
        error_message = 'Topology parsing error! (If you are sure everything is in order' \
                        ' try converting to another format)'
    except getattr(__builtins__,'FileNotFoundError', IOError):
        error_message = 'Topology file missing!'
    number_frames_total = 0

    if error_message == '':
        for trajectory in trajectories:
            try:
                system = MDAnalysis.Universe(topology, format=extension)
                number_frames_trajectory = len(system.trajectory)
                number_frames_total += number_frames_trajectory
            except ValueError:
                try:
                    system = MDAnalysis.Universe(topology, trajectory, format=extension)
                    number_frames_trajectory = len(system.trajectory)
                    number_frames_total += number_frames_trajectory
                except getattr(__builtins__, 'FileNotFoundError', IOError):
                    error_message = 'One or more files missing!'
                except (KeyError, ValueError) as error:
                    error_message = 'Topology/trajectory error or filetype mismatch!'

    if error_message == '':
        if ncores > multiprocessing.cpu_count():
            error_message = 'Number of cores specified higher than available number of cores (%s)'\
                            % multiprocessing.cpu_count()
        elif ncores > (number_frames_total // stride) + 1:
                error_message = 'Number of cores higher than number of frames!'

    if error_message == '':
        try:
            system.select_atoms(mask1)
        except (TypeError, KeyError, MDAnalysis.exceptions.SelectionError) as error:
            error_message = 'Mask1 error!'
        if error_message == '':
            try:
                system.select_atoms(mask2)
            except (TypeError, KeyError, MDAnalysis.exceptions.SelectionError) as error:
                error_message = 'Mask2 error!'

    if error_message == '':
        if len(system.select_atoms(mask1)) == 0:
            error_message = 'Mask1 contains no atoms!'
    if error_message == '':
        if len(system.select_atoms(mask2)) == 0:
            error_message = 'Mask2 contains no atoms!'
    if error_message == '':
        if cutoff != '':
            try:
                cutoff = float(cutoff)
            except ValueError:
                error_message = 'Cutoff must be a number!'
    if error_message == '' and calculate_vdw:
        try:
            system.atoms[0].element == 'a'
        except:
            error_message = 'Van der Waals interactions not supported for this format!'

    return error_message


def check_plot_params(settings, filename='config_plot.ini'):
    """Check validity of the parameters in the config_plot.ini file.
        Parameters
        ----------
        settings: ConfigParser
            Object containing data located in config_plot.ini
        filename: str
            Filename of the config file containing plotting parameters
        """
    error_message = ''
    datafile = settings.get('files', 'datafile')
    startframe = load_from_plot_config('parameters', 'startframe', filename)
    endframe = load_from_plot_config('parameters', 'endframe', filename)
    starting_residue = load_from_plot_config('parameters', 'starting_residue', filename)
    end_residue = load_from_plot_config('parameters', 'end_residue', filename)
    range_of_values = load_from_plot_config('parameters', 'range_of_values', filename)

    # Check filename
    try:
        with open(datafile, 'r+') as results:
            last_line = tail(results)
    except getattr(__builtins__, 'FileNotFoundError', IOError):
        error_message = 'File not found!'

    # Check datafile entry validity
    if error_message == '':
        try:
            final_datafile_frame = int(last_line[0].split()[0])
            last_datafile_residue = int(last_line[0].split()[1])
        except (ValueError, KeyError) as error:
            error_message = 'Datafile error!'

    # Check if parameters are integers

    try:
        startframe = int(startframe)
    except ValueError:
        error_message = 'Starting frame must be an integer!'

    try:
        endframe = int(endframe)
    except ValueError:
        error_message = 'Final frame must be an integer!'

    try:
        starting_residue = int(starting_residue)
    except ValueError:
        error_message = 'Starting residue must be an integer!'

    try:
        end_residue = int(end_residue)
    except ValueError:
        error_message = 'Final residue must be an integer!'

    if error_message == '':
        if startframe > endframe:
            error_message = 'Starting frame number higher than final frame!'

    if error_message == '':
        if startframe > final_datafile_frame:
            error_message = 'Starting frame larger than total number of frames (%d)' % final_datafile_frame
        elif starting_residue > last_datafile_residue:
                error_message = 'Starting residue larger than total number of residues (%d)' % last_datafile_residue

    if error_message == '':
        try:
            range_of_values = float(range_of_values)
        except ValueError:
            error_message = 'Range of values must be a number!'
    return error_message


def load_from_config(section, parameter, filename='config.ini'):
    """Tries loading parameters from the config.ini file. Defaults to values defined in default_values if not found in
    config.ini.

    Parameters
    ----------
    section: str
        Section of the .ini file
    parameter: str
        Parameter as defined in the .ini file
    filename: str
        Name of the config file

    Returns
    -------
        A list containing either the parameter from the input file or a default value and a boolean telling whether
        or not a default value was used (important for logging).
    """

    default_values = ConfigParser()
    config = ConfigParser()
    default_values.read(pkg_resources.resource_filename('relictoolkit', 'config_defaults.dft'))

    config.read(filename)

    try:
        if config.get(section, parameter) == '':
            return [default_values[parameter], True]
        else:
            return [config.get(section, parameter), False]
    except:
        return [default_values.get(section, parameter), True]


def load_from_plot_config(section, parameter, filename='config_plot.ini'):
    """Tries loading parameters from the config_plot.ini file. Defaults to values defined in default_values if not
    found in config_load.ini.

    Parameters
    ----------
    section: str
        Section of the .ini file
    parameter: str
        Parameter as defined in the .ini file
    filename: str
        Name of the config file
    Returns
    -------
        A parameter read from the config file or a default value defined in default_values
    """
    config = ConfigParser()
    config.read(filename)
    datafile = config.get('files', 'datafile')
    with open(datafile, 'r+') as results:
        last_line = tail(results)

    default_values = {
        'interactive': 'True',
        'plot_type': 'time',
        'startframe': 0,
        'endframe': int(last_line[0].split()[0]),
        'starting_residue': 1,
        'end_residue': int(last_line[0].split()[1]),
        'dt': 1,
        'range_of_values': 0
    }

    try:
        if config.get(section, parameter) == '':
            return default_values[parameter]
        else:
            return config.get(section, parameter)
    except:
        return default_values[parameter]


def load_partial_traj(input_system, step, ncores, core):
    """Perform trajectory pre-processing: apply periodic boundary conditions and break into multiple files for parallel
        processing.

        Parameters
        ----------
        input_system : MDAnalysis.Universe
            A Universe object containing the input system data loaded from topology/trajectory files
        step: int
            Trajectory sampling interval
        ncores: int
            Number of CPU cores
        core: int
            current CPU core
        """

    frames_per_traj = len(input_system.trajectory) // ncores

    # Overhead is used to avoid loading multiple frames within a single step due to the trajectory fragmenting
    # or improper fragmentation of the trajectory
    overhead = frames_per_traj - (frames_per_traj // step) * step
    overhead = overhead * core
    overhead = overhead % step
    endframe = (core + 1) * frames_per_traj
    if overhead == 0:
        overhead = step

    elif overhead == step:
        endframe += 1

    startingframe = (core*frames_per_traj) + (step-overhead)

    if core == ncores-1:
        endframe = len(input_system.trajectory)

    with open('relic_logfile.log', 'a+') as logfile:
        logfile.write('Core %s assigned frames %s to %s\n' % (core, startingframe, endframe))

    partial_traj_info = {
        'traj': input_system.trajectory[startingframe:endframe:step],
        'startframe': startingframe,
        'endframe': endframe
    }

    return partial_traj_info


def electrostatic_energy(atom1, atom2, box_dimensions):
    """Calculate electrostatic interactios within a cutoff distance.

    Parameters
    ----------
    atom1: MDAnalysis.core.groups.Atom
        First atom
    atom2: MDAnalysis.core.groups.Atom
        Second atom
    box_dimensions: numpy.ndarray
        Dimensions of the system box (for PBC imaging)

    Returns
    -------
    eelec: float
        Electrostatic energy in kJ/mol.
    """

    indiel = float(load_from_config('parameters', 'indi')[0])
    interatom_distance = MDAnalysis.lib.distances.calc_bonds(atom1.position.reshape(-1, 3),
                                                             atom2.position.reshape(-1, 3), box=box_dimensions)[0]
    eelec = 1389 * (1/(4*pi*indiel))*(atom1.charge * atom2.charge) / interatom_distance  # kJ/mol
    return eelec  # converts to kJ/mol


def interdomain_interactions(domain1, domain2, frame_number):
    """Calculate interactions between domain1 and domain2 present in the given frame

    Parameters
    ----------
    domain1: MDAnalysis.core.groups.AtomGroup
        Atoms comprising the first domain
    domain2: MDAnalysis.core.groups.AtomGroup
        Atoms comprising the second domain
    frame_number: int
        Frame number, numbering consistent with the original trajectory (for output)

    Returns
    -------
    domain_interactions: list
        List contaning frame number, residue number, electrostatic, van der Waals and total interaction energy. This
        list is subsequently printed to output.
    """

    possible_neighbors = AtomNeighborSearch(domain2.select_atoms('not resname WAT', updating=True))
    domain_interactions = list()
    calculate_vdw = load_from_config('parameters', 'vdw')[0] == 'True'
    for residue in domain1.residues:
        if residue.resname != 'WAT' and residue.resname != 'SOL':
            neighbor_atoms = possible_neighbors.search(residue.atoms,
                                                       float(load_from_config('parameters', 'cutoff')[0]), level='A')
            eelec = 0
            evdw = 0
            for residue_atom in residue.atoms:
                for neighbor_atom in neighbor_atoms:
                    eelec += electrostatic_energy(residue_atom, neighbor_atom, domain1.dimensions)
                    if calculate_vdw:
                        evdw += vdw_energy(residue_atom, neighbor_atom, domain1.dimensions)
            domain_interactions.append('{:<10d} {:<10d} {:>10.5f} {:>10.5f} {:>10.5f}'.format(
                frame_number, residue.resid, eelec, evdw, eelec+evdw))
    return domain_interactions


def process_frame(domain1, domain2, output, frame_number):
    """Calculate interactions for a single frame, sort by residue and print them to output.

    Parameters
    ----------
    domain1: MDAnalysis.core.groups.AtomGroup
        Atoms comprising the first domain
    domain2: MDAnalysis.core.groups.AtomGroup
        Atoms comprising the second domain
    output: TextIO
        Filename of the output file
    frame_number: int
        Frame number, numbering consistent with the original trajectory (for output)
    """

    frame_unsorted = list()
    interactions_domain1 = interdomain_interactions(domain1, domain2, frame_number)
    interactions_domain2 = interdomain_interactions(domain2, domain1, frame_number)
    frame_unsorted += interactions_domain1
    frame_unsorted += interactions_domain2
    for entry in sorted(frame_unsorted, key=lambda row: int(row.split()[1])):
        print(entry, file=output)


def read_uff_parameters():
    """Read LJ parameters from the UFF database. Based on Rappe et al., 1992.

    Returns
    -------
    parameters: dict
        Dictionary containing LJ parameters for every atom.
    """

    parameters = {}
    with open(pkg_resources.resource_filename('relictoolkit', 'uff.parm')) as uff_parameters:
        for line in uff_parameters:
            (key, val) = (line.split()[0], line.split()[1:3])
            try:
                val = [float(val[i]) for i in range(0, len(val))]
                parameters[key] = val
            except ValueError:
                pass
    return parameters


def write_logfile_header(settings, logfile):
    """Print the calculation start time and given input to the header of the logfile.
        Parameters
        ----------
        settings: ConfigParser
            Object containing data located in config.ini
        logfile: FileIO
            .log file"""

    logfile.write('Calculation started at %s\n\n' % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    logfile.write('-' * 15 + 'Input overview' + '-' * 15 + '\n')
    logfile.write('{:<17s}{:<80}\n'.format('Topology:', settings.get('files', 'topology')))
    trajectory_list = settings.get('files', 'trajectories').split(', ')
    logfile.write('Trajectories:\n')
    for name in trajectory_list:
        logfile.write('{:<17s}{:<80s}\n'.format('', name))
    logfile.write('Data written to: %s\n' % load_from_config('files', 'output')[0])
    logfile.write('Mask1: %s\nMask2: %s\n' % (settings.get('parameters', 'mask1'), settings.get('parameters', 'mask2')))

    strideinfo = load_from_config('parameters', 'stride')
    if not strideinfo[1]:
        logfile.write('Stride: %s\n' % int(strideinfo[0]))
    else:
        logfile.write('Stride not specified; a default value of %s used instead\n'
                      % int(strideinfo[0]))

    indiinfo = load_from_config('parameters', 'indi')
    if not indiinfo[1]:
        logfile.write('Internal dielectric constant: %s\n' % indiinfo[0])
    else:
        logfile.write('Internal dielectric constant not specified; using a default value of %s\n' % indiinfo[0])

    vdwinfo = load_from_config('parameters', 'vdw')
    if not vdwinfo[1]:
        if vdwinfo[0]:
            logfile.write('Van der Waals interactions will not be calculated\n')
        else:
            logfile.write('Van der Waals interactions will be calculated\n')
    else:
        logfile.write('Treatment of van der Waals interactions not specified; will not calculate\n')

    procinfo = load_from_config('parameters', 'ncores')
    if not procinfo[1]:
        logfile.write('Using %s CPU cores...\n' % procinfo[0])
    else:
        logfile.write('Number of cores not specified; using %s core\n' % procinfo[0])

    logfile.write('-' * 44 + '\n\n')


def tail(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end
    Taken from:
    https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail/7047765#7047765

    Parameters
    ----------
    f : fileI/O
        file containing the output
    lines : integer
        number of lines to get
    _buffer : integer
        chunk size

    Returns
    -------
        Selected number of lines from the end of the file
    """
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        # Removed this line because it was redundant the while will catch
        # it, I left it for history
        # if len(lines_found) > lines:
        #    break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]


def main():
    pass


if __name__ == '__main__':
    main()

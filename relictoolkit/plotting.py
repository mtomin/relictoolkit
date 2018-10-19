from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import int
from builtins import dict
from future import standard_library
import plotly.offline
import plotly.graph_objs as go
import relictoolkit.utils as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
standard_library.install_aliases()


def read_datafile(data, startframe, endframe, starting_residue, end_residue, plottype, dt):
    """Parse the calculation output and return datapoints for plotting.

    Parameters
    ----------
    data: _io.TextIOWrapper
        File containing the calculation output.
    startframe: int
        Starting frame for plotting
    endframe: int
        End frame for plotting
    starting_residue: int
        Starting residue for plotting
    end_residue: int
        End residue for pltoting
    plottype: str
        Plot type - energy vs. time, energy vs. frame number or average energy for each residue
    dt: int
        Timestep used in the simulations
    """

    frame = 0
    x = list()
    y = list()
    z = list()
    step = 1
    for line in data:
        current_frame = int(line.split()[0])
        current_res = int(line.split()[1])

        if current_frame < startframe or (current_res < starting_residue or current_res > end_residue):
            continue

        if frame != current_frame:

            step = current_frame - frame
            frame += step

            x.append(None)
            y.append(None)
            z.append(None)

        if plottype == 'time':
            x.append(float(line.split()[0])/(1000*dt))
        elif plottype == 'frame_number':
            x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))
        z.append(float(line.split()[2]))

        if frame > endframe:
            break

    datapoints = {'x': x,
                  'y': y,
                  'z': z,
                  'step': step}

    return datapoints


def generate_layout(plottype):
    if plottype == 'time':
        xaxis_title = 'Time/ns'
    elif plottype == 'frame_number':
        xaxis_title = 'Frame'
    elif plottype == 'averages':
        xaxis_title = 'Average energy'
    if plottype != 'averages':
        layout = go.Layout(
                title='Residue interactions',
                scene=dict(
                    xaxis=dict(title=xaxis_title,
                               autorange=True),
                    yaxis=dict(title='Residue',
                               autorange=True),
                    zaxis=dict(title='E/kJmol<sup>-1</sup>',
                               autorange=True),
                    aspectratio=dict(x=1, y=1, z=1),
                ),
                font=dict(family='Arial', size=18),
            )

    elif plottype == 'averages':
        layout = go.Layout(
            title='Average residue energies',
            scene=dict(
                xaxis=dict(title='Residue',
                           autorange=True),
                yaxis=dict(title='E<sub>avg</sub>/kJmol<sup>-1</sup>',
                           autorange=True),
            ),
            font=dict(family='Arial', size=18),
        )

    return layout


def calculate_averages(filename, datapoints, startframe, endframe, starting_residue, end_residue):
    with open(filename, 'r+') as results:
        last_line = u.tail(results)
    num_of_residues = int(last_line[0].split()[1])
    average_residue_energies = list()
    residues = list()

    for residue in range(starting_residue, end_residue + 1):
        avg_res_energy = 0
        for energy in datapoints['z'][residue - 1:len(datapoints['z']):num_of_residues + 1]:
            avg_res_energy += energy
        avg_res_energy = avg_res_energy * datapoints['step'] / (endframe - startframe + datapoints['step'])
        residues.append(residue)
        average_residue_energies.append(avg_res_energy)
    results = {
        'residues': residues,
        'residue_energies': average_residue_energies
    }
    return results


def generate_figure_data_plotly(config_filename):
    # Load parameters from config
    filename = u.load_from_plot_config('files', 'datafile', config_filename)
    plottype = u.load_from_plot_config('parameters', 'plot_type', config_filename)
    startframe = int(u.load_from_plot_config('parameters', 'startframe', config_filename))
    endframe = int(u.load_from_plot_config('parameters', 'endframe', config_filename))
    starting_residue = int(u.load_from_plot_config('parameters', 'starting_residue', config_filename))
    end_residue = int(u.load_from_plot_config('parameters', 'end_residue', config_filename))

    data = open(filename)

    traces = list()

    # Read timestep from file
    next(data)
    dt = int(data.readline().split()[1])

    datapoints = read_datafile(data, startframe, endframe, starting_residue, end_residue, plottype, dt)
    data.close()

    x = datapoints['x']
    y = datapoints['y']
    z = datapoints['z']

    traces.append(go.Scatter3d(
        mode='lines',
        line={'width': 5},
        x=x,
        y=y,
        z=z,
        connectgaps=False))

    if plottype != 'averages':
        fig = {
            'data': traces,
            'layout': generate_layout(plottype)
        }

    elif plottype == 'averages':
        averages_data = calculate_averages(filename, datapoints, startframe, endframe, starting_residue, end_residue)

        traces = [(go.Scatter(
            mode='lines',
            line={'width': 5},
            x=averages_data['residues'],
            y=averages_data['residue_energies'],
            connectgaps=False))]
        fig = {
            'data': traces,
            'layout': generate_layout(plottype)
        }

    return fig


def generate_figure_data_mplt(config_filename):
    # Load parameters from config
    filename = u.load_from_plot_config('files', 'datafile', config_filename)
    plottype = u.load_from_plot_config('parameters', 'plot_type', config_filename)
    startframe = int(u.load_from_plot_config('parameters', 'startframe', config_filename))
    endframe = int(u.load_from_plot_config('parameters', 'endframe', config_filename))
    starting_residue = int(u.load_from_plot_config('parameters', 'starting_residue', config_filename))
    end_residue = int(u.load_from_plot_config('parameters', 'end_residue', config_filename))

    data = open(filename)

    # Read timestep from file
    next(data)
    dt = int(data.readline().split()[1])

    datapoints = read_datafile(data, startframe, endframe, starting_residue, end_residue, plottype, dt)
    data.close()
    fig = plt.figure(facecolor='w')
    mpl.rcParams['savefig.dpi'] = 1000
    xcrds = list()
    ycrds = list()
    zcrds = list()
    for i,j,k in zip(datapoints['x'], datapoints['y'], datapoints['z']):
        print(i,j,k)
    if plottype == 'averages':
        averages_data = calculate_averages(filename, datapoints, startframe, endframe, starting_residue, end_residue)
        for x, y in zip(averages_data['residues'], averages_data['residue_energies']):
            xcrds.append(x)
            ycrds.append(y)

        ax = fig.add_subplot(111)
        ax.plot(xcrds, ycrds)
        ax.set_title('Average residue energies')
        ax.set_xlabel('Residue')
        ax.set_ylabel('E$_{avg}$/kJmol$^{-1}$')
        fig.tight_layout()
    else:
        ax = fig.add_subplot(111, projection='3d')
        for x, y, z in zip(datapoints['x'], datapoints['y'], datapoints['z']):

            if x is not None:
                xcrds.append(x)
                ycrds.append(y)
                zcrds.append(z)
            else:
                xnew, ynew = np.meshgrid(xcrds, ycrds)
                znew = np.zeros((len(ycrds), len(xcrds)))
                for i in range(len(zcrds)):
                    znew[i] = zcrds[i]
                xcrds = list()
                ycrds = list()
                zcrds = list()
                ax.plot_wireframe(xnew, ynew, znew, lw=0.25)

        ax.set_title('Residue interactions')
        if plottype == 'time':
            ax.set_xlabel('Time/ns')
        elif plottype == 'frame_number':
            ax.set_xlabel('Frame')
        ax.set_ylabel('Residue')
        ax.set_zlabel('E/kJmol$^{-1}$')
    ax.autoscale(enable=True, axis='both')
    return fig


def main(config_filename='config_plot.ini'):
    plottype = u.load_from_plot_config('parameters', 'plot_type', config_filename)
    interactive = u.load_from_plot_config('parameters', 'interactive', config_filename)
    if interactive == 'False':
        figure = generate_figure_data_mplt(config_filename)
        figure.show()
    else:
        figure = generate_figure_data_plotly(config_filename)
        if plottype != 'averages':
            plot_display_options = dict(toImageButtonOptions=dict(width=2400, height=2400, filename='relic_plot'))
            plotly.offline.plot(figure, auto_open=True, config=plot_display_options, filename='relic_plot.html')

        elif plottype == 'averages':
            plot_display_options = dict(toImageButtonOptions=dict(width=2400, height=2400, filename='relic_plot'))
            plotly.offline.plot(figure, auto_open=True, config=plot_display_options, filename='relic_plot.html')


if __name__ == '__main__':
    main()

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
standard_library.install_aliases()


def main(config_filename='config_plot.ini'):
    # Load parameters from config
    filename = u.load_from_plot_config('files', 'datafile', config_filename)
    plottype = u.load_from_plot_config('parameters', 'plot_type', config_filename)
    startframe = int(u.load_from_plot_config('parameters', 'startframe', config_filename))
    endframe = int(u.load_from_plot_config('parameters', 'endframe', config_filename))
    starting_residue = int(u.load_from_plot_config('parameters', 'starting_residue', config_filename))
    end_residue = int(u.load_from_plot_config('parameters', 'end_residue', config_filename))

    data = open(filename)
    frame = 0
    x = list()
    y = list()
    z = list()
    traces = list()

    # read timestep from file
    next(data)
    dt = int(data.readline().split()[1])

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

    traces.append(go.Scatter3d(
        mode='lines',
        line={'width': 5},
        x=x,
        y=y,
        z=z,
        connectgaps=False))

    if plottype == 'time':
        xaxis_title = 'Time/ns'
    elif plottype == 'frame_number':
        xaxis_title = 'Frame'
    elif plottype == 'averages':
        xaxis_title = 'Average energy'
    if plottype != 'averages':
        fig = {
            'data': traces,
            'layout': go.Layout(
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
        }
        plot_display_options = dict(toImageButtonOptions=dict(width=2400, height=2400, filename='relic_plot'))
        plotly.offline.plot(fig, auto_open=True, config=plot_display_options, filename='relic_plot.html')

    elif plottype == 'averages':
        with open(filename, 'r+') as results:
            last_line = u.tail(results)
        num_of_residues = int(last_line[0].split()[1])
        average_residue_energies = list()
        residues = list()

        for residue in range(starting_residue, end_residue + 1):
            avg_res_energy = 0
            for energy in traces[0]['z'][residue - 1:len(traces[0]['z']):num_of_residues + 1]:
                avg_res_energy += energy
            avg_res_energy = avg_res_energy * step / (endframe - startframe + step)
            residues.append(residue)
            average_residue_energies.append(avg_res_energy)

        traces=[(go.Scatter(
            mode='lines',
            line={'width': 5},
            x=residues,
            y=average_residue_energies,
            connectgaps=False))]
        fig = {
            'data': traces,
            'layout': go.Layout(
                title='Average residue energies',
                scene=dict(
                    xaxis=dict(title='Residue',
                               autorange=True),
                    yaxis=dict(title='E<sub>avg</sub>/kJmol<sup>-1</sup>',
                               autorange=True),
                ),
                font=dict(family='Arial', size=18),
            )
        }
        plot_display_options = dict(toImageButtonOptions=dict(width=2400, height=2400, filename='relic_plot'))
        plotly.offline.plot(fig, auto_open=True, config=plot_display_options, filename='relic_plot.html')


if __name__ == '__main__':
    main()

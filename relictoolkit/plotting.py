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
    timeorframe = int(u.load_from_plot_config('parameters', 'timeorframe', config_filename))
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

        if timeorframe == 0:
            x.append(float(line.split()[0])/(1000*dt))
        else:
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

    if timeorframe == 0:
        xaxis_title = 'Time/ns'
    else:
        xaxis_title = 'Frame'
    fig = {
        'data': traces,
        'layout': go.Layout(
            title='Residue interactions',
            margin=dict(l=0, r=0, b=0, t=0),
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


if __name__ == '__main__':
    main()

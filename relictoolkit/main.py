from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import object
from future import standard_library
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter.messagebox import showerror, showinfo
from configparser import ConfigParser
import re
import relictoolkit.calculation as calculation
import relictoolkit.plotting as plotting
from relictoolkit.utils import check_params
from relictoolkit.utils import check_plot_params
standard_library.install_aliases()


class TotalGui(object):
    """
    Container class for the main app.
    """

    class BrowseButton(tk.Button):
        """Class for data acquired through file browsing dialog: trajectory, topology, output.
        """

        def __init__(self, labeltext, tab):
            """
            Parameters
            ----------
            labeltext : str
                Text used for the button label
            tab : ttk.Frame
                Denotes which window tab contains the button
            """
            tk.Button.__init__(self, master=tab, text="Browse...", anchor='center', width=15)
            self.grid(pady=10, column=1)
            self.filename = tk.StringVar()
            self.label = tk.Label(tab, text=labeltext, width=30, anchor='center')
            self.label.grid(column=0)
            self.content_label = tk.Label(master=tab, textvariable=self.filename, anchor='w',
                                          justify='left', wraplength=500, width=65, relief="sunken")
            self.content_label.grid(column=2, padx=10)

    class DropDown(tk.OptionMenu):
        """Class for the drop-down menu (currently only used for specifying the file type).
        """

        def __init__(self, labeltext, choices, tab):
            """
            Parameters
            ----------
                labeltext : str
                    Text used for the button label
                choices : list of str
                    Possible options in the drop-down menu
                tab : ttk.Frame
                    Denotes which window tab contains the button
            """
            self.choices = choices
            self.currentval = tk.StringVar(tab)
            tk.OptionMenu.__init__(self, tab, self.currentval, *self.choices)
            self.label = tk.Label(tab, text=labeltext)
            self.currentval.set(choices[0])

    class EntryBox(tk.Entry):
        """Class for data entered through tk.Entry: masks, indi, stride, ncores.
        """

        def __init__(self, labeltext, tab):
            """
            Parameters
            ----------
            labeltext : str
                Text used for the button label
             tab : ttk.Frame
                    Denotes which window tab contains the button
            """
            tk.Entry.__init__(self, tab, width=15)
            self.label = tk.Label(tab, text=labeltext+':')
            self.label.grid(column=0)
            self.grid(column=1, pady=5)

    class CheckBox(tk.Checkbutton):
        """Class for data entered through tk.Checkbutton: calculating van der Waals interactions
        """

        def __init__(self, labeltext, tab):
            """
            Parameters
            ----------
            labeltext : str
                Text used for the button label
            """

            self.box_state = tk.BooleanVar()
            self.label = tk.Label(tab, text=labeltext)
            tk.Checkbutton.__init__(self, tab, variable=self.box_state)

    class RunButton(tk.Button):
        """Class for Run and Quit buttons
        """

        def __init__(self, labeltext, tab):
            """
            Parameters
            ----------
            labeltext : str
                Text used for the button label
             tab : ttk.Frame
                    Denotes which window tab contains the button
            """

            tk.Button.__init__(self, master=tab, text=labeltext)
            self.grid(pady=20, row=11)

    class Footer(tk.Label):
        """
        Class containing the signature footer.
        """

        def __init__(self, master):
            """
            Parameters
            ----------
            master : ttk.Frame
                tab containing the footer
            """

            tk.Label.__init__(self, master, height=1, bg='lightgrey', text='\N{COPYRIGHT SIGN} Marko Tomin', anchor='w')
            self.grid_configure(row=1, sticky='nsew')

    class Calculate (ttk.Frame):
        """The Calculate tab.
        """
        def __init__(self, master):
            """
            Parameters
            ----------
            master : tkinter.Tk()
                Master window of the Calculate frame
            """

            ttk.Frame.__init__(self, master)

            # Create browse buttons, labels and content labels for selecting trajectory, topology and output file
            self.trajectory = TotalGui.BrowseButton('Trajectories', self)
            self.topology = TotalGui.BrowseButton('Topology / Standalone format', self)

            self.output = TotalGui.BrowseButton('Set output path', self)
            self.topology.label.grid_configure(row=0)
            self.trajectory.label.grid_configure(row=1)
            self.output.label.grid_configure(row=2)

            # Assign the command to the topology browse button
            self.topology.content_label.grid_configure(row=0)
            self.topology.configure(command=lambda: self.topology.filename.set(browse_topology()))
            self.topology.grid_configure(row=0)

            # Assign the command to the trajectory browse button
            self.trajectory.content_label.grid_configure(row=1)
            self.trajectory.configure(command=lambda: self.trajectory.filename.set(browse_trajectory()))
            self.trajectory.grid_configure(row=1)

            # Assign the command to the output browse button
            self.output.content_label.grid_configure(row=2)
            self.output.configure(command=lambda: set_output(self.output))
            self.output.grid_configure(row=2)

            # Create format selection
            possible_filetypes = ['Amber binary trajectory (.ncdf)', 'AMBER coordinates (.inpcrd)',
                                  'Amber trajectory (.mdcrd)', 'CHARMM/XPLOR/PSF (.psf)', 'DESRES (.dms)',
                                  'GAMESS (.gms)', 'Gromacs trr (.trr)', 'Gromacs xtc (.xtc)', 'GROMOS96 (.gro)',
                                  'GSD (.gsd)', 'HOOMD (.xml)', 'LAMMPS (.data)']

            self.extension = TotalGui.DropDown('Filetype: ', possible_filetypes, self)
            self.extension.label.grid(row=3, column=1)
            self.extension.grid(row=3, column=2, sticky='w', padx=10)

            # Create entry boxes for masks, dielectric constant, stride, etc.
            self.mask1 = TotalGui.EntryBox('Mask1', self)
            self.mask2 = TotalGui.EntryBox('Mask2', self)
            self.cutoff = TotalGui.EntryBox('Distance cutoff / \N{ANGSTROM SIGN}', self)
            self.indi = TotalGui.EntryBox('Internal dielectric constant', self)
            self.stride = TotalGui.EntryBox('Stride', self)
            self.dt = TotalGui.EntryBox('Timestep / fs', self)
            self.processors = TotalGui.EntryBox('Number of cores', self)

            self.vdw = TotalGui.CheckBox('Calculate van der Waals interactions', self)

            # Chapter titled Options
            tk.Label(master=self, text='Options').grid(row=4, column=0, pady=10)

            # LAYOUT

            self.mask1.label.grid_configure(row=5)
            self.mask2.label.grid_configure(row=6)
            self.cutoff.label.grid_configure(row=7)
            self.indi.label.grid_configure(row=8)
            self.stride.label.grid_configure(row=9)
            self.dt.label.grid_configure(row=10)
            self.processors.label.grid_configure(row=5, column=2, sticky='w', padx=10)
            self.vdw.label.grid_configure(row=6, column=2, sticky='w', padx=10)

            self.mask1.grid_configure(row=5)
            self.mask2.grid_configure(row=6)
            self.cutoff.grid_configure(row=7)
            self.indi.grid_configure(row=8)
            self.stride.grid_configure(row=9)
            self.processors.grid_configure(row=5, column=2, sticky='w', padx=150)
            self.dt.grid_configure(row=10)

            self.vdw.grid_configure(row=6, column=2, sticky='w', padx=250)

            # Create Run and Quit buttons for the Calculate tab
            run_calc = TotalGui.RunButton('Run', self)
            run_calc.configure(command=lambda: run_calculation(self))
            run_calc.grid_configure(column=1)
            quit_program = TotalGui.RunButton('Quit', self)
            quit_program.configure(command=self.master.quit)
            quit_program.grid_configure(column=2)

    class Plot(ttk.Frame):
        """
        The Plot tab.
        """
        def __init__(self, master):
            """
            Parameters
            ----------
            master : tkinter.Tk()
                Master window of the Calculate frame
            """

            ttk.Frame.__init__(self, master)

            # Create and set up the data browsing button
            self.data = TotalGui.BrowseButton('Select datafile', self)
            self.data.label.grid_configure(row=0)
            self.data.grid_configure(row=0)
            self.data.content_label.grid_configure(row=0)
            self.data.filename.trace('w', self.generate_settings_pane)
            self.data.configure(text='Browse...', command=lambda: self.data.filename.set(browse_input()))

            self.plottype = tk.StringVar()
            self.interactive = tk.StringVar()
            self.startingframe = TotalGui.EntryBox('Starting frame (default 0)', self)
            self.endframe = TotalGui.EntryBox('Final frame (default last)', self)
            self.starting_residue = TotalGui.EntryBox('Starting residue (default 1)', self)
            self.end_residue = TotalGui.EntryBox('Final residue (default last)', self)

            # Hide settings until the data file is selected
            self.startingframe.grid_remove()
            self.endframe.grid_remove()
            self.starting_residue.grid_remove()
            self.end_residue.grid_remove()
            self.startingframe.label.grid_remove()
            self.endframe.label.grid_remove()
            self.starting_residue.label.grid_remove()
            self.end_residue.label.grid_remove()

        def generate_settings_pane(self, *args):
            """
            Spawns the plotting settings pane after the data file has been selected.

            Parameters
            ----------
            Arguments returned by trace
            """

            if self.data.filename.get() != 'None' and self.data.filename.get() != '':

                # Interactive
                interactive_label = tk.Label(self, text='Select plot drawing method:')
                interactive_label.grid(row=1, columnspan=3, sticky='n')
                matplotlib_radiobutton = ttk.Radiobutton(self, text='Non-interactive (matplotlib)',
                                                         variable=self.interactive,
                                                         value='False')
                matplotlib_radiobutton.grid_configure(row=2, columnspan=3, sticky='w', padx=250, pady=(10, 10))
                plotly_radiobutton = ttk.Radiobutton(self, text='Interactive (Plotly)',
                                                     variable=self.interactive,
                                                     value='True')
                plotly_radiobutton.grid_configure(row=2, columnspan=3, sticky='e', padx=250, pady=(10, 10))
                self.interactive.set('False')

                # Plot type layout
                plot_type_label = tk.Label(self, text='Plot type:')
                plot_type_label.grid(row=3, columnspan=3, sticky='n')

                time_radiobutton = ttk.Radiobutton(self,
                                                   text='Interactions vs. time',
                                                   variable=self.plottype,
                                                   value='time')
                time_radiobutton.grid_configure(row=4, columnspan=3, sticky='w', padx=100, pady=(5, 20))

                frame_radiobutton = ttk.Radiobutton(self,
                                                    text='Interactions vs. frame number',
                                                    variable=self.plottype,
                                                    value='frame_number')
                frame_radiobutton.grid_configure(row=4, columnspan=3, sticky='n', pady=(5, 20))

                averages_radiobutton = ttk.Radiobutton(self,
                                                       text='Average energy per residue',
                                                       variable=self.plottype,
                                                       value='averages')
                averages_radiobutton.grid_configure(row=4, columnspan=3, sticky='e', padx=(0, 100), pady=(5, 20))
                self.plottype.set('time')

                # Plot settings layout
                self.startingframe.grid_configure(row=5, column=1)
                self.startingframe.label.grid_configure(row=5, column=0)

                self.endframe.grid_configure(row=6, column=1)
                self.endframe.label.grid_configure(row=6, column=0)

                self.starting_residue.grid_configure(row=7, column=1)
                self.starting_residue.label.grid_configure(row=7, column=0)

                self.end_residue.grid_configure(row=8, column=1)
                self.end_residue.label.grid_configure(row=8, column=0)

                # Create and set up Run and Quit buttons
                run_plot = TotalGui.RunButton('Run', self)
                run_plot.configure(command=lambda: run_plotting(self))
                run_plot.grid_configure(column=1)
                quit_program = TotalGui.RunButton('Quit', self)
                quit_program.configure(command=self.master.quit)
                quit_program.grid_configure(column=2)

    def __init__(self, master):
        """
        Parameters
        ----------
        master : tkinter.Tk()
            Main window
        """
        self.master = master
        self.master.title('RELIC toolkit')
        self.master.resizable(False, False)
        self.Footer(master)

        # Create tabs Calculate and Plot
        notebook = ttk.Notebook(self.master)
        calculate = self.Calculate(self.master)
        plot = self.Plot(self.master)
        notebook.add(calculate, text='Calculate')
        notebook.add(plot, text='Plot')
        notebook.grid(row=0)


def browse_trajectory():
    """Allow user to select multiple trajectory files and update the trajectory label with the absolute path.

    Returns
    ----------
    filenames_formatted : list of str
        A formatted list of trajectory filenames
    """

    filelist = filedialog.askopenfilenames(filetypes=[('all files', '.*')], title='Open trajectory files')
    filenames_formatted = []
    for filename in filelist:
        filenames_formatted.append(filename)
    filenames_formatted = '\n'.join(filenames_formatted)
    return filenames_formatted


def browse_topology():
    """Allow user to select a single topology file and update the topology label with the absolute path

    Returns:
        filename : str
            Filename of the topology file
    """

    filename = filedialog.askopenfilenames(filetypes=[('all files', '.*')], title='Open topology file')

    # Check if a file has been selected
    if len(filename) > 0:
        return filename[0]


def set_output(out):
    """Prompt user to select the destination for the output file and update the output label with the absolute path.

    Parameters
    ----------
    out : TotalGui.BrowseButton
        Contains the filename variable that is updated in the method
    """

    filename = filedialog.asksaveasfilename(filetypes=[('all files', '.*')], title='Set output file destination')
    if len(filename) > 0:
        out.filename.set(filename)


def browse_input():
    """Prompt user to select the file containing data to be plotted and update the output label with the absolute path.
    """

    filename = filedialog.askopenfilenames(filetypes=[('all files', '.*')], title='Open input file')

    # Check if a file has been selected
    if len(filename) > 0:
        return filename[0]


def generate_config(gui, config):
    """After the input has been entered, create a config.ini to be used for calculation

    Parameters
    ----------
    gui : tkinter.Tk()
        Interface
    config : configparser.Configparser
        Container for data found in config.ini
    """

    # Get extension from the drop-down menu entry
    ext = re.search(r'\(.*\)', gui.extension.currentval.get()).group()

    config['files'] = {
        'topology': gui.topology.filename.get(),
        'trajectories': gui.trajectory.filename.get().replace('\n', ', '),
        'output': gui.output.filename.get(),
        'filetype': ext[2:len(ext)-1]
    }
    # Set values in the 'parameters' section of the config file
    config['parameters'] = {
        'mask1': gui.mask1.get(),
        'mask2': gui.mask2.get(),
        'indi': gui.indi.get(),
        'vdw': gui.vdw.box_state.get(),
        'stride': gui.stride.get(),
        'cutoff': gui.cutoff.get(),
        'ncores': gui.processors.get(),
        'dt': gui.dt.get()
    }
    with open('config.ini', 'w+') as configfile:
        config.write(configfile)


def generate_plot_config(gui, config):
    """
    Take data from GUI and write it to config_plot.ini.

    Parameters
    ----------
    gui : Plot(ttk.Frame)
        Frame containing the data to be written to config
    config : ConfigParser
        Configparser object containing the data taken from GUI
    """

    config['files'] = {
        'datafile': gui.data.filename.get()
    }
    config['parameters'] = {
        'interactive': gui.interactive.get(),
        'plot_type': gui.plottype.get(),
        'startframe': gui.startingframe.get(),
        'endframe': gui.endframe.get(),
        'starting_residue': gui.starting_residue.get(),
        'end_residue': gui.end_residue.get()
    }

    with open('config_plot.ini', 'w+') as configfile:
        config.write(configfile)


def run_calculation(gui):
    """Make the config.ini file for calculating electrostatic interactions, check for input errors and initialize
    the electrostatic interaction calculation.

    Parameters
    ----------
    gui : Calculate(ttk.Frame)
        Frame containing the data to be written to config
    """

    config = ConfigParser()
    generate_config(gui, config)
    error_message = check_params(config)
    if error_message == '':
        calculation.main()
        showinfo('Done', 'Calculation completed successfully!')
    else:
        showerror('Error', error_message)


def run_plotting(gui):
    """Call plotting.py to plot the given dataset

    Parameters
    ----------
    gui : ttk.Frame
        Frame containing the data to be written to config
    """
    plot_config = ConfigParser()
    generate_plot_config(gui, plot_config)
    error_message = check_plot_params(plot_config)
    if error_message == '':
        plotting.main()
    else:
        showerror('Error', error_message)


def main():
    root = tk.Tk()
    relic_gui = TotalGui(root)
    root.mainloop()


if __name__ == "__main__":
    main()

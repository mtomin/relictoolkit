from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object
from future import standard_library
import relictoolkit.main as r
import tkinter as tk
from tkinter import ttk
from configparser import ConfigParser
import os
standard_library.install_aliases()

os.environ['DISPLAY'] = '0.0'


class TestPlot(ttk.Frame):
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
        self.data = r.TotalGui.BrowseButton('Select datafile', self)
        self.data.label.grid_configure(row=0)
        self.data.grid_configure(row=0)
        self.data.content_label.grid_configure(row=0)
        self.data.configure(text='Browse...')

        self.plottype = tk.StringVar()
        self.startingframe = r.TotalGui.EntryBox('Starting frame (default 0)', self)
        self.endframe = r.TotalGui.EntryBox('Final frame (default last)', self)
        self.starting_residue = r.TotalGui.EntryBox('Starting residue (default 1)', self)
        self.end_residue = r.TotalGui.EntryBox('Final residue (default last)', self)

        # Hide settings until the data file is selected
        self.startingframe.grid_remove()
        self.endframe.grid_remove()
        self.starting_residue.grid_remove()
        self.end_residue.grid_remove()
        self.startingframe.label.grid_remove()
        self.endframe.label.grid_remove()
        self.starting_residue.label.grid_remove()
        self.end_residue.label.grid_remove()


class TestGui(object):
    def __init__(self):
        testtab = ttk.Frame()

        # For testing generate_config
        self.topology = TestGui.TestGuivarfiles('topology')
        self.trajectory = TestGui.TestGuivarfiles('trajectory')
        self.output = TestGui.TestGuivarfiles('output')
        self.mask1 = TestGui.TestGuivarparameters('mask1')
        self.mask2 = TestGui.TestGuivarparameters('mask2')
        self.indi = TestGui.TestGuivarparameters('indi')
        self.stride = TestGui.TestGuivarparameters('stride')
        self.cutoff = TestGui.TestGuivarparameters('cutoff')
        self.processors = TestGui.TestGuivarparameters('processors')
        self.dt = TestGui.TestGuivarparameters('dt')
        self.extension = r.TotalGui.DropDown('extensions', ['Amber binary trajectory (.ncdf)'], testtab)
        self.vdw = TestGui.TestGuivarparameters('vdw')

        # For testing generate_plot_config
        self.data = TestGui.TestGuivarfiles('datafile')
        self.plottype = TestGui.TestGuivarparameters('plottype')
        self.startingframe = TestGui.TestGuivarparameters('startframe')
        self.endframe = TestGui.TestGuivarparameters('endframe')
        self.starting_residue = TestGui.TestGuivarparameters('starting_residue')
        self.end_residue = TestGui.TestGuivarparameters('end_residue')

    class TestGuivarfiles(tk.StringVar):
        def __init__(self, value):
            tk.StringVar.__init__(self)
            self.filename = tk.StringVar()
            self.filename.set(value)

    class TestGuivarparameters(tk.StringVar):
        def __init__(self, value):
            tk.StringVar.__init__(self, value=value)
            self.box_state = tk.StringVar()
            self.box_state.set(value)


def test_totalgui_browsebutton():
    testtab = ttk.Frame()
    testbutton = r.TotalGui.BrowseButton('testbutton', testtab)
    assert testbutton.grid_info()['pady'] == 10
    assert testbutton.label.cget('text') == 'testbutton'


def test_totalgui_dropdown():
    testtab = ttk.Frame()
    choices = ['choice1', 'choice2', 'choice3']
    testmenu = r.TotalGui.DropDown('testmenu', choices, testtab)
    assert testmenu.currentval.get() == 'choice1'
    assert testmenu.label.cget('text') == 'testmenu'


def test_totalgui_entrybox():
    testtab = ttk.Frame()
    testentrybox = r.TotalGui.EntryBox('testentrybox', testtab)
    assert testentrybox.label.cget('text') == 'testentrybox:'
    assert testentrybox.grid_info()['pady'] == 5


def test_totalgui_checkbox():
    testtab = ttk.Frame()
    testcheckbox = r.TotalGui.CheckBox('testcheckbox', testtab)
    assert testcheckbox.label.cget('text') == 'testcheckbox'


def test_totalgui_runbutton():
    testtab = ttk.Frame()
    testrunbutton = r.TotalGui.RunButton('testrunbutton', testtab)
    assert testrunbutton.cget('text') == 'testrunbutton'
    assert testrunbutton.grid_info()['pady'] == 20


def test_totalgui_footer():
    testtab = ttk.Frame()
    testfooter = r.TotalGui.Footer(testtab)
    assert testfooter.cget('bg') == 'lightgrey'


def test_totalgui_calculate():
    testroot = tk.Tk()
    testcalculate = r.TotalGui.Calculate(testroot)
    assert testcalculate.trajectory.label.cget('text') == 'Trajectories'


def test_totalgui_plot():
    testroot = tk.Tk()
    testplot = r.TotalGui.Plot(testroot)
    assert testplot.data.label.cget('text') == 'Select datafile'


def test_totalgui_plot_generate_settings_pane():
    testroot = tk.Tk()
    testplot = TestPlot(testroot)
    r.TotalGui.Plot.generate_settings_pane(testplot)
    assert testplot.plottype.get() == ''
    testplot.data.filename.set('randomfilename')
    assert testplot.plottype.get() == ''
    r.TotalGui.Plot.generate_settings_pane(testplot)
    assert testplot.plottype.get() == 'time'


def test_generate_config():
    testgui = TestGui()
    config = ConfigParser()
    r.generate_config(testgui, config)
    assert config.get('files', 'filetype') == 'ncdf'
    assert config.get('files', 'topology') == 'topology'
    assert config.get('parameters', 'ncores') == 'processors'
    os.remove('config.ini')


def test_generate_plot_config():
    testgui = TestGui()
    config = ConfigParser()
    r.generate_plot_config(testgui, config)
    assert config.get('files', 'datafile') == 'datafile'
    assert config.get('parameters', 'plot_type') == 'plottype'
    assert config.get('parameters', 'end_residue') == 'end_residue'
    os.remove('config_plot.ini')

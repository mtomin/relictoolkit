# Ignore warnings due to multiple warnings regarding numpy in dependencies etc.
import sys
import os
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from setuptools.command.install import install as InstallCommand

class PyTest(TestCommand):
    """
    Unit test wrapper for the PyTest, including coverage repport
    """
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--color=yes', 'tests/'] #['--cov spectraplotpy tests/']
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

with open('README.rst') as f:
    readme = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='patkica',
    version='0.9',
    setup_requires=['numpy',
		    'matplotlib==1.5.2',
		],
    packages=find_packages(),
    zip_safe=True,
    scripts=[],
    include_package_data=True,
    install_requires=[
	required
    ],
    
    entry_points = {
        'console_scripts': [
            'relictoolkit = relictoolkit.main:main',
	    'relictoolkit_calculate = relictoolkit.calculation:main',
	    'relictoolkit_plot = relictoolkit.plotting:main'
                 ],              
          },

    package_data={
        '': 'relictoolkit',
	},

    # Project uses pytest for the tests
    tests_require=[
        'pytest',
        'pytest-cov',
        'mock'
    ],
    
    cmdclass={
        'test': PyTest,
    },

    url="www.irb.hr",
    license="CC BY 4.0",
    author="Marko Tomin",
    author_email="marko.tomin@irb.hr",
    description="Residue ELectrostatic Interaction Calculator",
    long_description=readme,
)

#ipywidgets==7.4.2

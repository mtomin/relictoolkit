from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from future import standard_library
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    """
    Unit test wrapper for the PyTest, including coverage repport
    """
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--color=yes', 'tests/'] # ['--cov spectraplotpy tests/']
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


with open('README.md') as f:
    readme = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='relictoolkit',
    version='0.9',
    setup_requires=['setuptools>=40.0.0',
                    'numpy',
                    'matplotlib==1.5.2'],
    packages=find_packages(),
    zip_safe=True,
    scripts=[],
    include_package_data=True,
    install_requires=[required],
    
    entry_points={
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

    #url='RELIC github, https://github.com/mtomin/relictoolkit',
    license="CC BY 4.0",
    author="Marko Tomin",
    author_email="marko.tomin@irb.hr",
    description="Residue ELectrostatic Interaction Calculator",
    long_description=readme,
    long_description_content_type="text/markdown",
)

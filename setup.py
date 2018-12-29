#!/usr/bin/env python

""" Setup for da_controller """

from glob import glob
from setuptools import find_packages, setup
from os.path import basename, dirname, realpath

def main():
    return setup(
        author='Forecasting Research group',
        author_email='research@metservice.com',
        data_files=[('var_blending/etc', glob('etc/*'))],
        description='to be added',
        maintainer='Forecasting Research group',
        maintainer_email='research@metservice.com',
        # Use name of the directory as name
        name=basename(dirname(realpath(__file__))),
        scripts=['scripts/run_var_blending.py',
                 'scripts/run_var_blending_err.py'],
        packages=find_packages(),
        zip_safe=False
    )


if __name__ == '__main__':
    main()

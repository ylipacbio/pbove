from setuptools import setup, Extension, find_packages
from distutils.extension import Extension


setup(
    name = 'pbove',
    version='0.1.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license='LICENSE.txt',
    scripts = ['pbove/pbove_eval.py',
               'pbove/pbove_compare_runs.py',
               'pbove/pbove_compare_overlap.py',
               'pbove/pbove_plot_runs.py',
               'pbove/pbove_filter_subreads.py',
               'pbove/pbove_print_overlaps.py',
               'pbove/pbove_reseq.py',
               'pbove/pbove_preassembly.py'],
    package_dir={'pbove':'pbove'},
    package_data={'pbove':['R/pbove_compare_runs.R']},
    packages = find_packages(),
    entry_points={'console_scripts':
                  ['pbove=pbove.pbove_main:main']},
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.6.3',
        'pbalign >= 0.1.0']
    )

# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

# General imports
import os
import sys
import time as t
import json
import argparse
from subprocess import run

# Custom imports
sys.path.append('src')
from plot_data import plot_all

# ----------------------------------------------------------------------------------------------------------------------
# Program Execution
# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Runs Einstein-aether collapse simulation.')
    parser.add_argument('-o', type=str, default='results', metavar='output',
                        help='Path to simulation output and graphs. Default value is results.')
    parser.add_argument('-p', type=str, default='param/params.json', metavar='params',
                        help='Path to .json file containing simulation parameters. Default value is param/params.json.')
    args = parser.parse_args()

    # Initialize paths for saving data and plots
    save_data = f'{args.o}/data/'
    save_plots = f'{args.o}/plots/'

    # Creating paths if they don't already exist
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    
    if not os.path.exists(save_data):
        os.mkdir(save_data)

    if not os.path.exists(save_plots):
        os.mkdir(save_plots)

    # Get simulation parameters from .json file
    with open(args.p) as f:
        params = json.load(f)
    f.close()

    amp1 = params['amp1']
    amp2 = params['amp2']
    r0 = params['r0']
    sigma = params['sigma']
    rmax = params['rmax']
    ntime = params['ntime']
    time = params['time']
    tmax = params['tmax']
    c1 = params['c1']
    c2 = params['c2']
    c3 = params['c3']
    c4 = params['c4']

    # Run collapse simulation
    print('Running simulation...')
    t0 = t.time()

    compile_sim = ['gfortran', 'src/aether.f90', '-o', 'src/aether.exe']
    run(compile_sim)

    run_sim = ['src/aether.exe', save_data, f'{amp1}', f'{amp2}', f'{r0}', f'{sigma}', f'{rmax}', f'{ntime}', f'{time}',
               f'{tmax}', f'{c1}', f'{c2}', f'{c3}', f'{c4}']
    run(run_sim)

    tf = t.time()
    print(f'Simulation completed in {int((tf - t0)//60)} min {round((tf - t0) % 60, 1)} sec.')
    print('\n')

    # Plot simulation results
    plot_all(save_data, save_plots)

# ----------------------------------------------------------------------------------------------------------------------
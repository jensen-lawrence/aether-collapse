# Spherically Symmetric Collapse in Einstein-Aether Theory #

## Introduction ##
The code in this repository simulates spherically symmetric collapse in Einstein-aether theory with an additional scalar matter field. Specifically, Crank-Nicholson iteration is used to time-evolve Eqs. (C2, C3, C5-C10) in [2].

## Software Used ##
The simulation code is written in Fortran, using post-Fortran 90 conventions. A Fortran compiler (such as [gfortran](https://gcc.gnu.org/wiki/GFortran)) is required to compile the code. The data processing and plot creation code is written in Python, version 3.7 or higher.

In addition, the Python portion of the code makes use of the following libraries:
- os
- sys
- time
- json
- argparse
- subprocess
- numpy
- matplotlib
- colour

## Usage ##
To run this program, navigate to the `aether-collapse` directory on the command line, i.e., `cd path/to/aether-collapse`. Then, run the command
```
python main.py -o output -p params
```
Note that on computers with Python 2 and Python 3 installed, it may be necessary to use `python3` instead of `python` in the previous command.

The first flag, `-o`, is mandatory. Its argument, `output`, is a string. This is the directory in which the simulation data and plots will be saved. Its default value is `results`. If an alternate value is provided as the argument, this will be used as the output directory instead.

The second flag, `-p`, is mandatory. Its argument, `params`, is a string. This is the path to the `.json` file containing the simulation parameters. Its default value is `param/params.json`. If an alternate value is provided as the argument, this will be used as the parameters file instead.

## References ##
[1] https://arxiv.org/abs/gr-qc/0703093

[2] https://arxiv.org/abs/1608.06970

# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from colour import Color

# ----------------------------------------------------------------------------------------------------------------------
# Comparison Plot Function
# ----------------------------------------------------------------------------------------------------------------------

def compare_results(results, labels, output):
    """
    Docstring goes here
    """
    # Output data file names
    filenames = ['alphafin', 'arfin', 'arfin2', 'asqfin', 'cnstrfin', 'cnstrfin2', 'curvfin', 'grrfin', 'kfin', 'kfin2',
                 'massfin', 'pfin', 'psifin', 'qfin2', 'radfin', 'tnnfin', 'trapfin', 'trapfin2', 'trapmin', 'trapmin2']

    # Plot labels
    xlabels = [r'$r$', r'$r$', r'$R$', r'$R$', r'$r$', r'$r$', r'$\log{z}$', r'$R$', r'$r$', r'$R$', r'$R$', r'$r$',
               r'$r$', r'$R$', r'$r$', r'$r$', r'$r$', r'$r$', r'$t$', r'$t$']
    ylabels = [r'$\alpha$', r'$a_r$', r'$a_r$', r'$A^2$', r'$C_1$', r'$C_2$', r'$\log(\mathrm{Curv})$', r'$g_{rr}$',
               r'$K$', r'$K$', r'$M$', r'$P$', r'$\psi$', r'$Q$', r'$R$', r'$T_{nn}$', 'Trapped Surface',
               'Anti-Trapped Surface', 'Minimum Trapped Surface Value', 'Minimum Anti-Trapped Surface Value']

    # Plot colours
    min_colour = Color('#00FF7F')
    max_colour = Color('#8A2BE2')
    colour_range = list(min_colour.range_to(max_colour, len(results)))

    # Generate plots
    for i in range(len(filenames)):
        filename = filenames[i]
        xlabel = xlabels[i]
        ylabel = ylabels[i]

        fig, ax = plt.subplots(figsize=(16,9))
        ax.grid()

        for j in range(len(results)):
            dir = results[j]
            colour = colour_range[j]
            label = labels[j]

            data = np.loadtxt(dir + f'/{filename}.txt')
            x, y = data[:,0], data[:,1]
            ax.plot(x, y, color=str(colour), lw=2, label=label)

        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.legend(loc='best', fontsize=16)
        plt.savefig(f'{output}/{filename}_comp.png', dpi=300)
        plt.close()

# ----------------------------------------------------------------------------------------------------------------------
# Program Execution
# ----------------------------------------------------------------------------------------------------------------------

# Results to be plotted
results = [ 
    '../results/maxv0/c1=0.01_c2=0.01/data',
    '../results/maxv0/c1=0.01_c2=0.05/data',
    '../results/maxv0/c1=0.01_c2=0.10/data',
    '../results/maxv0/c1=0.01_c2=0.15/data',
    '../results/maxv0/c1=0.01_c2=0.20/data',
    '../results/maxv0/c1=0.01_c2=0.30/data',
    '../results/maxv0/c1=0.01_c2=0.40/data'
]

labels = [ 
    r'$c_1 = 0.01, \; c_2 = 0.01$',
    r'$c_1 = 0.01, \; c_2 = 0.05$',
    r'$c_1 = 0.01, \; c_2 = 0.10$',
    r'$c_1 = 0.01, \; c_2 = 0.15$',
    r'$c_1 = 0.01, \; c_2 = 0.20$',
    r'$c_1 = 0.01, \; c_2 = 0.30$',
    r'$c_1 = 0.01, \; c_2 = 0.40$',
]

output = '/home/jlawrence/aether-collapse/results/maxv0/comparison'

# compare_results(results, labels, output)

vvals = np.array([1.4038116913620973, 2.3629671849069869, 3.0938502827237087, 3.6157748419908629, 4.0215950216368270,
                  4.3515209011477189, 4.6274666056756226, 4.8629463473967034, 5.0669702344132972])

tvals = np.array([366.8, 619.7, 807, 940.9, 1056.6, 1153.5, 1209.1, 1270.8, 1328.4])

fig, ax = plt.subplots(figsize=(16,9))
ax.grid(True)
ax.plot(vvals, tvals, color='dodgerblue', lw=2)
ax.set_xlabel(r'$v_0$', fontsize=16)
ax.set_ylabel(r'$t$', fontsize=16)
plt.savefig('t_vs_v0.png', dpi=300)
plt.close()

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import animation

# ----------------------------------------------------------------------------------------------------------------------
# Create Static Plot
# ----------------------------------------------------------------------------------------------------------------------

def static_plot(data_file, xlabel, ylabel, output, dpi=300):
    """
    Docstring goes here
    """
    # Load data
    data = np.loadtxt(data_file)
    savelabel = data_file.split('/')[-1].split('.')[0]

    # Get x and y values
    xvals = data[:,0]
    yvals = data[:,1]

    # Get plotting limits for x and y values
    xdiff = 0.05*(np.amax(xvals) - np.amin(xvals))
    ydiff = 0.05*(np.amax(yvals) - np.amin(yvals))
    xlims = (np.amin(xvals) - xdiff, np.amax(xvals) + xdiff)
    ylims = (np.amin(yvals) - ydiff, np.amax(yvals) + ydiff)

    # Initialize plot
    fig, ax = plt.subplots(figsize=(16,9))
    ax.grid()
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    # Plot and save results
    ax.plot(xvals, yvals, color='dodgerblue', lw=2)
    plt.savefig(f'{output}/{savelabel}.png', dpi=dpi)
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------
# Create Animated Plot
# ----------------------------------------------------------------------------------------------------------------------

def animated_plot(data_file, xlabel, ylabel, output, fps=15):
    """
    Docstring goes here
    """
    # Load data
    data = np.loadtxt(data_file)
    savelabel = data_file.split('/')[-1].split('.')[0]

    # Get time values
    all_t = data[:,0]
    tvals = np.unique(all_t)
    
    # Get x and y values
    xvals, yvals = [], []
    for t in tvals:
        data_at_t = data[data[:,0] == t]
        x = data_at_t[:,1]
        y = data_at_t[:,2]
        xvals.append(x)
        yvals.append(y)
    xvals, yvals = np.array(xvals), np.array(yvals) 

    # Get plotting limits for x and y values
    xdiff = 0.05*(np.amax(xvals) - np.amin(xvals))
    ydiff = 0.05*(np.amax(yvals) - np.amin(yvals))
    xlims = (np.amin(xvals) - xdiff, np.amax(xvals) + xdiff)
    ylims = (np.amin(yvals) - ydiff, np.amax(yvals) + ydiff)

    # Initialize plot
    fig, ax = plt.subplots(figsize=(16,9))
    ax.grid()
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    # Initialize plot elements to be animated
    line, = ax.plot([], [], color='dodgerblue', lw=2)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=16)

    # Initialization function for animator
    def init():
        line.set_data([], [])
        return line, 

    # Animation function for animator
    def anim(i):
        x = xvals[i]
        y = yvals[i]
        line.set_data(x, y)
        time_text.set_text(f'{round(tvals[i], 3)}')
        return line, time_text

    # Animate and save results
    anim_len = len(tvals)
    anim_result = animation.FuncAnimation(fig, anim, init_func=init, frames=anim_len, interval=500, blit=True)
    writergif = animation.PillowWriter(fps=fps)
    anim_result.save(f'{output}/{savelabel}.gif', writer=writergif)
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------
# Plot All Variables
# ----------------------------------------------------------------------------------------------------------------------

def plot_all(data_in, plots_out):
    """
    Docstring goes here
    """
    # Generate static plots
    print('Generating static plots...')
    static_plot(f'{data_in}/alphafin.txt', r'$r$', r'$\alpha$', plots_out)
    static_plot(f'{data_in}/arfin.txt', r'$r$', r'$a_r$', plots_out)
    static_plot(f'{data_in}/arfin2.txt', r'$R$', r'$a_r$', plots_out)
    static_plot(f'{data_in}/asqfin.txt', r'$R$', r'$A^2$', plots_out)
    static_plot(f'{data_in}/cnstrfin.txt', r'$r$', r'$C_1$', plots_out)
    static_plot(f'{data_in}/cnstrfin2.txt', r'$r$', r'$C_2$', plots_out)
    static_plot(f'{data_in}/curvfin.txt', r'$\log{R}$', r'$\log(\mathrm{Curv})$', plots_out)
    static_plot(f'{data_in}/grrfin.txt', r'$R$', r'$g_{rr}$', plots_out)
    static_plot(f'{data_in}/kfin.txt', r'$r$', r'$K$', plots_out)
    static_plot(f'{data_in}/kfin2.txt', r'$R$', r'$K$', plots_out)
    static_plot(f'{data_in}/massfin.txt', r'$R$', r'$M$', plots_out)
    static_plot(f'{data_in}/pfin.txt', r'$r$', r'$P$', plots_out)
    static_plot(f'{data_in}/psifin.txt', r'$r$', r'$\psi$', plots_out)
    static_plot(f'{data_in}/qfin2.txt', r'$R$', r'$Q$', plots_out)
    static_plot(f'{data_in}/radfin.txt', r'$r$', r'$R$', plots_out)
    static_plot(f'{data_in}/tnnfin.txt', r'$r$', r'$T_{nn}$', plots_out)
    static_plot(f'{data_in}/trapfin.txt', r'$r$', 'Trapped Surface', plots_out)
    static_plot(f'{data_in}/trapfin2.txt', r'$r$', 'Anti-Trapped Surface', plots_out)
    static_plot(f'{data_in}/trapmin.txt', r'$t$', 'Minimum Trapped Surface Value', plots_out)
    static_plot(f'{data_in}/trapmin2.txt', r'$t$', 'Minimum Anti-Trapped Surface Value', plots_out)

    # Generate animated plots
    print('Generating animated plots...')
    animated_plot(f'{data_in}/Psi.txt', r'$r$', r'$\psi$', plots_out)
    animated_plot(f'{data_in}/P.txt', r'$r$', r'$P$', plots_out)
    animated_plot(f'{data_in}/ar.txt', r'$r$', r'$a_r$', plots_out)
    animated_plot(f'{data_in}/K.txt', r'$r$', r'$K$', plots_out)
    animated_plot(f'{data_in}/alpha.txt', r'$r$', r'$\alpha$', plots_out)
    animated_plot(f'{data_in}/Tll.txt', r'$r$', r'$T_{ll}$', plots_out)
    animated_plot(f'{data_in}/Tnn.txt', r'$r$', r'$T_{nn}$', plots_out)
    print('Plotting complete.')

# ----------------------------------------------------------------------------------------------------------------------
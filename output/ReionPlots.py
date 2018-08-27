#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

from astropy import cosmology

from scipy import stats

import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import ObservationalData as Obs

from mpi4py import MPI

output_format = "png"

def plot_history(rank, comm, z_array_reion_allmodels, 
                 lookback_array_reion_allmodels, mass_frac_allmodels, cosmo,
                 t_bigbang, model_tags, output_dir, output_tag):

    # Mass fractions are effectively histograms (single number for each z).
    master_mass_frac= collective.collect_hist_across_tasks(rank, comm, 
                                                           mass_frac_allmodels)

    if rank != 0:
        return 

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_reion_allmodels)):

        ax1.plot(lookback_array_reion_allmodels[model_number], 
                 master_mass_frac[model_number], 
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

    ax1.set_xlabel(r"$\mathbf{Time \: since \: Big \: Bang \: [Myr]}$", 
                   fontsize = ps.global_labelsize)
    tick_locs = np.arange(200.0, 1000.0, 100.0)
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(ps.time_tickinterval))
    tick_labels = [r"$\mathbf{%d}$" % x for x in tick_locs]
    ax1.xaxis.set_major_locator(mtick.MultipleLocator(100))
    ax1.set_xticklabels(tick_labels, fontsize = ps.global_fontsize)

    ax1.set_ylabel(r'$\mathbf{\langle \chi_{HI} \rangle}$',
                   fontsize = ps.global_labelsize)

    tick_locs = np.arange(-0.2, 1.1, 0.2)
    ax1.set_yticklabels([r"$\mathbf{%.2f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))
    ax1.set_ylim([-0.05, 1.05])

    ax1.set_xlim(ps.time_xlim)

    ax2 = ps.add_time_z_axis(ax1, cosmo[0], t_bigbang[0]/1.0e3)

    for axis in ['top','bottom','left','right']: # Adjust axis thickness.
        ax1.spines[axis].set_linewidth(ps.global_axiswidth)

    ax1.tick_params(which = 'both', direction='in', width = ps.global_tickwidth)
    ax1.tick_params(which = 'major', length = ps.global_ticklength)
    ax1.tick_params(which = 'minor', length = ps.global_ticklength-2)

    leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)

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
import ObservationalData as obs

from mpi4py import MPI


def plot_nion(z_array_full_allmodels, lookback_array_full_allmodels,
              cosmology_allmodels, t_bigbang_allmodels, nion_allmodels,
              model_tags, output_dir, output_tag, output_format):
    """
    Plots the ionizing emissivity (using the galaxy data files).

    Parameters
    ----------

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    lookback_array_reion_allmodels : 2D nested list of floats. Dimensions
                                     identical to ``z_array_reion_allmodels``. 
        The lookback time at each snapshot for each model. Units are Myr.

    cosmology_allmodels : List of class ``astropy.cosmology``. Length is number
                          of models.
        ``astropy`` class containing the cosmology for each model.

    t_bigbang_allmodels : List of floats. Length is number of models.
        The lookback time to the Big Bang for each model. Units is Gyr.

    nion_allmodels : 2D nested list of floats. Dimensions equal to
                     ``z_array_reion_allmodels``.
        The number of ionizing photons (that escaped) for each model at each
        snapshot. Units is number per Mpc^3.

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String.
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String.
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_full_allmodels)):

        # Units for photons are 1.0e50 so move to log properly.
        log_sum = 50 + np.log10(nion_allmodels[model_number])

        ax1.plot(lookback_array_full_allmodels[model_number], 
                 log_sum, 
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

    ax1 = ps.plot_bouwens2015(cosmology_allmodels[0], t_bigbang_allmodels[0],
                              ax1)

    ax1.set_xlabel(r"$\mathbf{Time \: since \: Big \: Bang \: [Myr]}$", 
                   fontsize = ps.global_labelsize)
    tick_locs = np.arange(200.0, 1000.0, 100.0)
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(ps.time_tickinterval))
    tick_labels = [r"$\mathbf{%d}$" % x for x in tick_locs]
    ax1.xaxis.set_major_locator(mtick.MultipleLocator(100))
    ax1.set_xticklabels(tick_labels, fontsize = ps.global_fontsize)

    ax1.set_ylabel(r'$\mathbf{log_{10} \sum f_{esc}\dot{N}_\gamma \: [s^{-1} Mpc^{-3}]}$',
                   fontsize = ps.global_labelsize)
    tick_locs = np.arange(48.0, 51.25, 0.5)
    ax1.set_yticklabels([r"$\mathbf{%.1f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)
    ax1.set_ylim([48.4, 51.25])
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.25))

    ax1.set_xlim(ps.time_xlim)

    ax2 = ps.add_time_z_axis(ax1, cosmology_allmodels[0],
                             t_bigbang_allmodels[0]/1.0e3)

    ax1.text(350, 50.0, r"$\mathbf{68\%}$", horizontalalignment='center', 
             verticalalignment = 'center', fontsize = ps.global_labelsize-4)
    ax1.text(350, 50.5, r"$\mathbf{95\%}$", horizontalalignment='center',
             verticalalignment = 'center', fontsize = ps.global_labelsize-4)

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


def plot_mstar_fesc(mstar_bins, mstar_bin_width, 
                    z_array_full_allmodels, mean_allmodels, std_allmodels,
                    N_allmodels, model_tags, output_dir, output_tag,
                    output_format, plot_models_at_snaps=None,
                    plot_snaps_for_models=None):
    """
    Plots the escape fraction as a function of stellar mass. 

    Parameters
    ----------

    mstar_bins : List of floats
        Stellar mass bins that the data is binned on. Units are Msun.

    mstar_bin_width : Float
        The bin separation between the stellar mass bins. 

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    mean_allmodels, std_allmodels, N_allmodels : 3D nested lists of floats.
                                                 Outer length is number of
                                                 models, next is number of
                                                 snapshots in the model and
                                                 final is the number of
                                                 ``mstar_bins``. 
        The mean and standard deviation for the escape fraction in each stellar
        mass bin.  Also the number of data points in each bin. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    plot_models_at_snaps : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots each model at the specified snapshots. That is,
        each panel will be for one model at the specified snapshots.

    plot_snaps_for_models : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots all models at a single, specified snapshot. That
        is, each panel will be for all models at one specified snapshot. 

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    # Our plotting area will be a square so find out the max NxN deimsnion.
    if plot_models_at_snaps: 
        nrows = int(np.ceil(np.sqrt(len(plot_models_at_snaps))))
    else:
        nrows = int(np.ceil(np.sqrt(len(plot_snaps_for_models)))) 

    fig, ax = plt.subplots(nrows=nrows, ncols=nrows, 
                           sharex='col', sharey='row', figsize=(16,6))

    row_ax = -1
    for count, model_number in enumerate(range(len(z_array_full_allmodels))):
        if count % nrows == 0:
            row_ax += 1
        col_ax = count % nrows
        for snap_count, snap in enumerate(plot_snaps_for_models[model_number]):

            z_label = r"$\mathbf{z = " + \
                        str(int(round(z_array_full_allmodels[model_number][snap],0))) + \
                       "}$"                

            mean = mean_allmodels[model_number][snap]
            w_low = np.where(master_N[model_number][snap] < 4)[0]
            mean[w_low] = np.nan

            ax[row_ax, col_ax]. plot(mstar_bins[:-1] + mstar_bin_width*0.5,
                                     mean, 
                                     color = ps.colors[snap_count],
                                     ls = ps.linestyles[0],
                                     label = z_label) 

        ax[row_ax, col_ax].text(0.05, 0.65, model_tags[model_number],
                                transform = ax[row_ax, col_ax].transAxes,
                                fontsize = ps.global_fontsize) 

        for axis in ['top','bottom','left','right']: # Adjust axis thickness.
            ax[row_ax, col_ax].spines[axis].set_linewidth(ps.global_axiswidth)

        ax[row_ax, col_ax].tick_params(which = 'both', direction='in', width =
                                       ps.global_tickwidth)
        ax[row_ax, col_ax].tick_params(which = 'major', length = ps.global_ticklength)
        ax[row_ax, col_ax].tick_params(which = 'minor', length = ps.global_ticklength-2)

    # Set variables for every column.
    tick_locs = np.arange(4.0, 11.0)
    for ax_count in range(nrows):
        ax[nrows-1, ax_count].set_xlabel(r'$\mathbf{log_{10} \: M_{*} \:[M_{\odot}]}$', 
                                         size = ps.global_fontsize)
        ax[nrows-1, ax_count].set_xlim([4.8, 10.2])
        ax[nrows-1, ax_count].xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax[nrows-1, ax_count].xaxis.set_major_locator(mtick.MultipleLocator(1.0))
        ax[nrows-1, ax_count].set_xticklabels([r"$\mathbf{%d}$" % x for x in tick_locs], 
                                              fontsize = ps.global_fontsize)
    
    labels = ax[1,0].xaxis.get_ticklabels()
    locs = ax[1,0].xaxis.get_ticklocs()
    for label, loc in zip(labels, locs):
        print("{0} {1}".format(label, loc)) 

    # Set variables for every row.
    tick_locs = np.arange(-0.10, 0.80, 0.10)
    for ax_count in range(nrows):
        ax[ax_count, 0].set_ylabel(r'$\mathbf{\langle f_{esc}\rangle_{M_*}}$', 
                                   size = ps.global_labelsize)

        ax[ax_count, 0].set_ylim([-0.02, 0.72])
        ax[ax_count, 0].yaxis.set_minor_locator(mtick.MultipleLocator(0.05))       
        ax[ax_count, 0].yaxis.set_major_locator(mtick.MultipleLocator(0.1))       
        ax[ax_count, 0].set_yticklabels([r"$\mathbf{%.2f}$" % x for x in tick_locs], 
                                         fontsize = ps.global_fontsize)

    leg = ax[0,0].legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig)


def plot_SMF(mstar_bins, mstar_bin_width, z_array_full_allmodels,
             cosmology_allmodels, SMF, model_tags, output_dir, output_tag,
             output_format, plot_z=[6.0, 7.0, 8.0]):
    """
    Plots the stellar mass function. That is, the number count of galaxies
    binned on stellar mass. 

    Parameters
    ----------

    mstar_bins : List of floats
        Stellar mass bins that the data is binned on. Units are Msun.

    mstar_bin_width : Float
        The bin separation between the stellar mass bins. 

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    cosmology_allmodels : List of class ``astropy.cosmology``. Length is number
                          of models.
        ``astropy`` class containing the cosmology for each model.

    SMF : 3D nested lists of floats. Outer length is number of models, next is
          number of snapshots in the model and final is the number of
          ``mstar_bins``.
        The stellar mass function at each snapshot for each model.  That is,
        the number of galaxies within each stellar mass bin
        (given by ``mstar_bins``).

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    plot_z : List of floats, optional
        The redshift at which we wish to plot the stellar mass function.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig, ax = plt.subplots(nrows=1, ncols=len(plot_z), 
                           sharex='col', sharey='row', figsize=(16,6))

    for model_number in range(len(z_array_full_allmodels)):
        # Here we find the snapshots closest to the redshifts we want to plot.
        plot_snaps = []
        for val in plot_z:
            idx = (np.abs(z_array_full_allmodels[model_number] - val)).argmin()
            plot_snaps.append(idx)

        # Then plot only those snapshots.
        for count, snap in enumerate(plot_snaps):
            label = model_tags[model_number]
            ax[count].plot(mstar_bins[:-1] + mstar_bin_width*0.5,
                           SMF[model_number][snap],
                           color = ps.colors[model_number],
                           ls = ps.linestyles[model_number],
                           label = label)

            if model_number == 0:
            
                tick_locs = np.arange(6.0, 12.0)
                ax[count].set_xticklabels([r"$\mathbf{%d}$" % x for x in tick_locs],
                                          fontsize = ps.global_fontsize)
                ax[count].set_xlim([6.8, 10.3])
                ax[count].set_xlabel(r'$\mathbf{log_{10} \: M_{*} \:[M_{\odot}]}$', 
                                     fontsize = ps.global_labelsize)
                ax[count].xaxis.set_minor_locator(plt.MultipleLocator(0.25))

                ax[count] = ps.adjust_axis(ax[count], ps.global_axiswidth,
                                           ps.global_tickwidth,
                                           ps.global_major_ticklength,
                                           ps.global_minor_ticklength) 

    # Since y-axis is shared, only need to do this once.
    ax[0].set_yscale('log', nonposy='clip')
    ax[0].set_yticklabels([r"$\mathbf{10^{-5}}$",r"$\mathbf{10^{-5}}$",r"$\mathbf{10^{-4}}$", r"$\mathbf{10^{-3}}$",
                           r"$\mathbf{10^{-2}}$",r"$\mathbf{10^{-1}}$"]) 
    ax[0].set_ylim([1e-5, 3e-1])
    #ax[0].set_ylabel(r'\mathbf{$\log_{10} \Phi\ [\mathrm{Mpc}^{-3}\: \mathrm{dex}^{-1}]}$', 
    ax[0].set_ylabel(r'$\mathbf{log_{10} \: \Phi\ [Mpc^{-3}\: dex^{-1}]}$', 
                     fontsize = ps.global_labelsize) 

    # Now lets overplot the Observational Data. 
    obs.Get_Data_SMF()

    caps = 5
    ewidth = 1.5

    ps.Plot_SMF_z6(ax[0], cosmology_allmodels[0].H(0).value/100.0, errorwidth=ewidth, capsize=caps) 
    ps.Plot_SMF_z7(ax[1], cosmology_allmodels[0].H(0).value/100.0, errorwidth=ewidth, capsize=caps) 
    ps.Plot_SMF_z8(ax[2], cosmology_allmodels[0].H(0).value/100.0, errorwidth=ewidth, capsize=caps) 
    
    ####

    delta_fontsize = 0
    ax[0].text(0.7, 0.9, r"$\mathbf{z = 6}$", transform = ax[0].transAxes, fontsize = ps.global_fontsize - delta_fontsize)
    ax[1].text(0.7, 0.9, r"$\mathbf{z = 7}$", transform = ax[1].transAxes, fontsize = ps.global_fontsize - delta_fontsize)
    ax[2].text(0.7, 0.9, r"$\mathbf{z = 8}$", transform = ax[2].transAxes, fontsize = ps.global_fontsize - delta_fontsize)

    leg = ax[0].legend(loc='lower left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize - 2)

    plt.tight_layout()

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig)


def plot_mstar_fej(mstar_bins, mstar_bin_width, z_array_full_allmodels,
                   mean_allmodels, std_allmodels, N_allmodels, model_tags,
                   output_dir, output_tag, output_format,
                   plot_models_at_snaps=None, plot_snaps_for_models=None):
    """
    Plots the fraction of baryons in the ejected reservoir as a function of
    stellar mass.

    Parameters
    ----------

    mstar_bins : List of floats
        Stellar mass bins that the data is binned on. Units are Msun.

    mstar_bin_width : Float
        The bin separation between the stellar mass bins. 

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.


    mean_allmodels, std_allmodels, N_allmodels : 3D nested lists of floats.
                                                 Outer length is number of
                                                 models, next is number of
                                                 snapshots in the model and
                                                 final is the number of
                                                 ``mstar_bins``. 
        The mean and standard deviation for the ejected fraction in each stellar
        mass bin.  Also the number of data points in each bin. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    plot_models_at_snaps : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots each model at the specified snapshots. That is,
        each panel will be for one model at the specified snapshots.

    plot_snaps_for_models : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots all models at a single, specified snapshot. That
        is, each panel will be for all models at one specified snapshot. 

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig, ax, nrows = plot_2D_line(mstar_bins, mstar_bin_width,
                                  z_array_full_allmodels, 
                                  mean_allmodels, std_allmodels, N_allmodels,
                                  model_tags, plot_models_at_snaps,
                                  plot_snaps_for_models)

    # Set variables for every column.
    tick_locs = np.arange(4.0, 11.0)
    for ax_count in range(nrows):
        ax[nrows-1, ax_count].set_xlabel(r'$\mathbf{log_{10} \: M_{*} \:[M_{\odot}]}$', 
                                         size = ps.global_fontsize)
        ax[nrows-1, ax_count].set_xlim([4.8, 10.2])
        ax[nrows-1, ax_count].xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax[nrows-1, ax_count].xaxis.set_major_locator(mtick.MultipleLocator(1.0))
        ax[nrows-1, ax_count].set_xticklabels([r"$\mathbf{%d}$" % x for x in tick_locs], 
                                              fontsize = ps.global_fontsize)
    
    labels = ax[1,0].xaxis.get_ticklabels()
    locs = ax[1,0].xaxis.get_ticklocs()
    for label, loc in zip(labels, locs):
        print("{0} {1}".format(label, loc)) 

    # Set variables for every row.
    tick_locs = np.arange(-0.10, 0.80, 0.10)
    for ax_count in range(nrows):
        ax[ax_count, 0].set_ylabel(r'$\mathbf{\langle f_{ej}\rangle_{M_*}}$', 
                                   size = ps.global_labelsize)

        ax[ax_count, 0].set_ylim([-0.02, 1.0])
        ax[ax_count, 0].yaxis.set_minor_locator(mtick.MultipleLocator(0.05))       
        ax[ax_count, 0].yaxis.set_major_locator(mtick.MultipleLocator(0.1))       
        #ax[ax_count, 0].set_yticklabels([r"$\mathbf{%.2f}$" % x for x in tick_locs], 
        #                                 fontsize = ps.global_fontsize)

    leg = ax[0,0].legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig)


def plot_mstar_SFR(mstar_bins, mstar_bin_width, z_array_full_allmodels,
                   mean_allmodels, std_allmodels, N_allmodels, model_tags,
                   output_dir, output_tag, output_format,
                   plot_models_at_snaps=None, plot_snaps_for_models=None):
    """
    Plots the star formation rate as a function of stellar mass.

    Parameters
    ----------

    mstar_bins : List of floats
        Stellar mass bins that the data is binned on. Units are Msun.

    mstar_bin_width : Float
        The bin separation between the stellar mass bins. 

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.


    mean_allmodels, std_allmodels, N_allmodels : 3D nested lists of floats.
                                                 Outer length is number of
                                                 models, next is number of
                                                 snapshots in the model and
                                                 final is the number of
                                                 ``mstar_bins``. 
        The mean and standard deviation for the SFR in each stellar mass bin.
        Also the number of data points in each bin. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    plot_models_at_snaps : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots each model at the specified snapshots. That is,
        each panel will be for one model at the specified snapshots.

    plot_snaps_for_models : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots all models at a single, specified snapshot. That
        is, each panel will be for all models at one specified snapshot. 

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig, ax, nrows = plot_2D_line(mstar_bins, mstar_bin_width,
                                  z_array_full_allmodels, mean_allmodels,
                                  std_allmodels, N_allmodels, model_tags,
                                  plot_models_at_snaps, plot_snaps_for_models)

    delta_fontsize = 10

    # Set variables for every column.
    tick_locs = np.arange(4.0, 11.0)
    for ax_count in range(nrows):
        ax[nrows-1, ax_count].set_xlabel(r'$\mathbf{log_{10} \: M_{*} \:[M_{\odot}]}$', 
                                         size = ps.global_fontsize)
        ax[nrows-1, ax_count].set_xlim([4.8, 10.2])
        ax[nrows-1, ax_count].xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax[nrows-1, ax_count].xaxis.set_major_locator(mtick.MultipleLocator(1.0))
        ax[nrows-1, ax_count].set_xticklabels([r"$\mathbf{%d}$" % x for x in tick_locs], 
                                              fontsize = ps.global_fontsize)
    
    labels = ax[1,0].xaxis.get_ticklabels()
    locs = ax[1,0].xaxis.get_ticklocs()
    for label, loc in zip(labels, locs):
        print("{0} {1}".format(label, loc)) 

    # Set variables for every row.
    tick_locs = np.arange(-0.10, 0.80, 0.10)
    for ax_count in range(nrows):
        ax[ax_count, 0].set_ylabel(r'$\mathbf{\langle SFR\rangle_{M_*} [M_\odot \: yr^{-1}]}$', 
                                   size = ps.global_labelsize-delta_fontsize)
        ax[ax_count, 0].set_yscale('log')
        #ax[ax_count, 0].set_ylim([-0.02, 1.0])
        #ax[ax_count, 0].yaxis.set_minor_locator(mtick.MultipleLocator(0.05))       
        #ax[ax_count, 0].yaxis.set_major_locator(mtick.MultipleLocator(0.1))       
        #ax[ax_count, 0].set_yticklabels([r"$\mathbf{%.2f}$" % x for x in tick_locs], 
        #                                 fontsize = ps.global_fontsize)

    leg = ax[0,0].legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig)


def plot_2D_line(bins, bin_width, binning_array_allmodels,
                 mean_array_allmodels, std_array_allmodels, N_array_allmodels,
                 model_tags, plot_models_at_snaps, plot_snaps_for_models):
    """
    Takes 2D binned arrays and creates a square plot of them. The number of
    rows and columns is equal to the number of models. 

    Parameters
    ----------

    bins : List of floats
        Bins that the data is binned on.

    bin_width : Float
        The bin separation between the bins. 

    mean_allmodels, std_allmodels, N_allmodels : 3D nested lists of floats.
                                                 Outer length is number of
                                                 models, next is number of
                                                 snapshots in the model and
                                                 final is the number of
                                                 ``bins``. 
        The mean and standard deviation for data in each bin. Also the
        number of data points in each bin. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    plot_models_at_snaps : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots each model at the specified snapshots. That is,
        each panel will be for one model at the specified snapshots.

    plot_snaps_for_models : 2D nested list of integers.  Outer length is number
                           of models, optional
        If not ``None``, plots all models at a single, specified snapshot. That
        is, each panel will be for all models at one specified snapshot. 

    Returns
    ---------

    fig : ``matplotlib`` Figure
        The created figure.

    ax : If number models is 1, a single ``matplotlib`` axes. Otherwise, 2D
         nested list of ``matplotlib`` axes with square dimension equal to the
         number of models.
        The axes for each panel in the figure.

    nrows : Integer
        The number of rows (and columns) in the figure.
    """

    # Our plotting area will be a square so find out the max NxN deimsnion.
    if plot_models_at_snaps: 
        nrows = int(np.ceil(np.sqrt(len(plot_models_at_snaps))))
    else:
        nrows = int(np.ceil(np.sqrt(len(plot_snaps_for_models)))) 

    fig, ax = plt.subplots(nrows=nrows, ncols=nrows, 
                           sharex='col', sharey='row', figsize=(16,6))

    row_ax = -1
    for count, model_number in enumerate(range(len(binning_array_allmodels))):
        if count % nrows == 0:
            row_ax += 1
        col_ax = count % nrows
        if nrows != 1:
            this_ax = ax[row_ax, col_ax]
        else:
            this_ax = ax

        for snap_count, snap in enumerate(plot_snaps_for_models[model_number]):

            z_label = r"$\mathbf{z = " + \
                        str(int(round(binning_array_allmodels[model_number][snap],0))) + \
                       "}$"

            mean = mean_array_allmodels[model_number][snap]
            w_low = np.where(N_array_allmodels[model_number][snap] < 4)[0]
            mean[w_low] = np.nan

            this_ax.plot(bins[:-1] + bin_width*0.5,
                         mean,
                         color = ps.colors[snap_count],
                         ls = ps.linestyles[0],
                         label = z_label)

        this_ax.text(0.05, 0.65, model_tags[model_number],
                     transform = this_ax.transAxes,
                     fontsize = ps.global_fontsize)

        this_ax = ps.adjust_axis(this_ax, ps.global_axiswidth,
                                 ps.global_tickwidth,
                                 ps.global_major_ticklength,
                                 ps.global_minor_ticklength)

    return fig, ax, nrows

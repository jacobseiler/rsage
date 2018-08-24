
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

def plot_zreion_hist(rank, comm, 
                     zreion_hist_allmodels, z_array_allmodels, model_tags,
                     output_dir, output_tag): 

    master_counts = collective.collect_hist_across_tasks(rank, comm, 
                                                         zreion_hist_allmodels)

    if rank != 0:
        return 

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_allmodels)):
        ax1.plot(z_array_allmodels[model_number], 
                 master_counts[model_number],
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

    ax1.set_yscale("log")

    ax1.set_xlabel(r"$\mathbf{z_{reion}}$",
                   fontsize = ps.global_labelsize)
    ax1.set_ylabel(r"$\mathbf{Count}$",
                   fontsize = ps.global_labelsize)

    leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


def plot_zreion_SFR(rank, comm, 
                    mean_zreion_SFR_allmodels, std_zreion_SFR_allmodels,
                    N_zreion_SFR_allmodels, zreion_bins_allmodels,
                    zreion_array_allmodels,
                    model_tags, output_dir, output_tag):   

    master_mean, master_std, master_N, _ = \
        collective.collect_across_tasks(rank, comm, \
                                        mean_zreion_SFR_allmodels, \
                                        std_zreion_SFR_allmodels, \
                                        N_zreion_SFR_allmodels, \
                                        zreion_bins_allmodels, \
                                        zreion_bins_allmodels, \
                                        binned=True)

    if rank != 0:
        return 

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)


    for model_number in range(len(zreion_bins_allmodels)):
        ax1.plot(zreion_array_allmodels[model_number], 
                 master_mean[model_number][10],
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

        z_label = r"$\mathbf{z = " + \
                    str(round(zreion_array_allmodels[model_number][10],2)) + \
                   "}$"                
        ax1.text(12.5, 1e-2, z_label, fontsize=ps.global_labelsize-4) 

    ax1.set_yscale("log")

    ax1.set_xlabel(r"$\mathbf{z}$",
                   fontsize = ps.global_labelsize)
    ax1.set_ylabel(r"$\mathbf{Mean SFR}$",
                   fontsize = ps.global_labelsize)

    leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


def plot_nion(rank, comm,
              z_array_full_allmodels, lookback_array_full_allmodels,
              sum_nion_allmodels, cosmo, t_bigbang, 
              model_tags, output_dir, output_tag):

    master_sum = collective.collect_hist_across_tasks(rank, comm, 
                                                      sum_nion_allmodels)

    if rank != 0:
        return 

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_full_allmodels)):

        # Units for photons are 1.0e50 so move to log properly.
        log_sum = 50 + np.log10(master_sum[model_number])

        ax1.plot(lookback_array_full_allmodels[model_number], 
                 log_sum, 
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])


    bouwens_z = np.arange(6,16) # Redshift range for the observations.
    bouwens_t = (t_bigbang[0] - cosmo[0].lookback_time(bouwens_z).value*1.0e3) # Corresponding values for what we will plot on the x-axis.

    bouwens_1sigma_lower = [50.81, 50.73, 50.60, 50.41, 50.21, 50.00, 49.80, 49.60, 49.39, 49.18] # 68% Confidence Intervals for the ionizing emissitivity from Bouwens 2015.
    bouwens_1sigma_upper = [51.04, 50.85, 50.71, 50.62, 50.56, 50.49, 50.43, 50.36, 50.29, 50.23]

    bouwens_2sigma_lower = [50.72, 50.69, 50.52, 50.27, 50.01, 49.75, 49.51, 49.24, 48.99, 48.74] # 95% CI.
    bouwens_2sigma_upper = [51.11, 50.90, 50.74, 50.69, 50.66, 50.64, 50.61, 50.59, 50.57, 50.55]
    
    ax1.fill_between(bouwens_t, bouwens_1sigma_lower, bouwens_1sigma_upper,
                     color = 'k', alpha = 0.7,
                     label = r"$\mathbf{Bouwens \: et \: al. \: (2015)}$")
    ax1.fill_between(bouwens_t, bouwens_2sigma_lower, bouwens_1sigma_lower,
                     color = 'k', hatch = '//', edgecolor = 'k', alpha = 0.4, 
                     facecolor = 'k', lw = 0.0)
    ax1.fill_between(bouwens_t, bouwens_1sigma_upper, bouwens_2sigma_upper,
                     color = 'k', hatch = '//', edgecolor = 'k', alpha = 0.4,
                     facecolor = 'k', lw = 0.0)

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

    ax2 = ps.add_time_z_axis(ax1, cosmo[0], t_bigbang[0]/1.0e3)

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


def plot_mstar_fesc(rank, comm, mstar_bins, mstar_bin_width, 
                    z_array_full_allmodels, mean_mstar_fesc_allmodels, 
                    std_mstar_fesc_allmodels, N_mstar_fesc_allmodels,
                    model_tags, output_dir, output_tag,
                    plot_models_at_snaps=None, plot_snaps_for_models=None):

    master_mean, master_std, master_N, _ = \
        collective.collect_across_tasks(rank, comm, \
                                        mean_mstar_fesc_allmodels, \
                                        std_mstar_fesc_allmodels, \
                                        N_mstar_fesc_allmodels, \
                                        z_array_full_allmodels, \
                                        z_array_full_allmodels, \
                                        binned=True)

    if rank != 0:
        return 

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

            mean = master_mean[model_number][snap]
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


def plot_SMF(rank, comm, mstar_bins, mstar_bin_width,
             z_array_full_allmodels, SMF, cosmology_allmodels, 
             model_tags, output_dir, output_tag,
             plot_z=[6.0, 7.0, 8.0]):

    master_counts = collective.collect_hist_across_tasks(rank, comm, 
                                                         SMF)

    if rank != 0:
        return 

    fig, ax = plt.subplots(nrows=1, ncols=len(plot_z), 
                           sharex='col', sharey='row', figsize=(16,6))

    for model_number in range(len(z_array_full_allmodels)):
        plot_snaps = []
        for val in plot_z:
            idx = (np.abs(z_array_full_allmodels[model_number] - val)).argmin()
            plot_snaps.append(idx)
        for count, snap in enumerate(plot_snaps):
            # If only one model is used, don't use a labels.
            if len(z_array_full_allmodels) == 1:
                label = ""
            else:
                label = model_tags[model_number]
            ax[count].plot(mstar_bins[:-1] + mstar_bin_width*0.5,
                             master_counts[model_number][snap],
                             color = ps.colors[model_number],
                             ls = ps.linestyles[model_number],
                             label = label)

            if model_number == 0:
            
                tick_locs = np.arange(6.0, 12.0)
                ax[count].set_xticklabels([r"$\mathbf{%d}$" % x for x in tick_locs],
                                          fontsize = ps.global_fontsize)
                ax[count].set_xlim([6.8, 10.3])                    
                ax[count].tick_params(which = 'both', direction='in', 
                                      width = ps.global_tickwidth)
                ax[count].tick_params(which = 'major',
                                      length = ps.global_ticklength)
                ax[count].tick_params(which = 'minor',
                                      length = ps.global_ticklength-2)
                ax[count].set_xlabel(r'$\mathbf{log_{10} \: M_{*} \:[M_{\odot}]}$', 
                                     fontsize = ps.global_labelsize)
                ax[count].xaxis.set_minor_locator(plt.MultipleLocator(0.25))
                #ax[count].set_xticks(np.arange(6.0, 12.0))
                
                for axis in ['top','bottom','left','right']: # Adjust axis thickness.
                    ax[count].spines[axis].set_linewidth(ps.global_axiswidth)

    # Since y-axis is shared, only need to do this once.
    ax[0].set_yscale('log', nonposy='clip')
    ax[0].set_yticklabels([r"$\mathbf{10^{-5}}$",r"$\mathbf{10^{-5}}$",r"$\mathbf{10^{-4}}$", r"$\mathbf{10^{-3}}$",
                           r"$\mathbf{10^{-2}}$",r"$\mathbf{10^{-1}}$"]) 
    ax[0].set_ylim([1e-5, 1e-1])
    #ax[0].set_ylabel(r'\mathbf{$\log_{10} \Phi\ [\mathrm{Mpc}^{-3}\: \mathrm{dex}^{-1}]}$', 
    ax[0].set_ylabel(r'$\mathbf{log_{10} \: \Phi\ [Mpc^{-3}\: dex^{-1}]}$', 
                     fontsize = ps.global_labelsize) 

    # Now lets overplot the Observational Data. 
    Obs.Get_Data_SMF()

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
                               
    #leg = ax[0,0].legend(loc=2, bbox_to_anchor = (0.2, -0.5), numpoints=1, labelspacing=0.1)
    leg = ax[0].legend(loc='lower left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize - 2)

    plt.tight_layout()

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig)

def plot_mstar_fej(rank, comm, mstar_bins, mstar_bin_width, 
                   z_array_full_allmodels, mean_mstar_fej_allmodels, 
                   std_mstar_fej_allmodels, N_mstar_fej_allmodels,
                   model_tags, output_dir, output_tag,
                   plot_models_at_snaps=None, plot_snaps_for_models=None):

    master_mean, master_std, master_N, _ = \
        collective.collect_across_tasks(rank, comm, \
                                        mean_mstar_fej_allmodels, \
                                        std_mstar_fej_allmodels, \
                                        N_mstar_fej_allmodels, \
                                        z_array_full_allmodels, \
                                        z_array_full_allmodels, \
                                        binned=True)

    if rank != 0:
        return 

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

            mean = master_mean[model_number][snap]
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


#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patheffects as PathEffects
import numpy as np

from astropy import cosmology

from scipy import stats

import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import ObservationalData as Obs

from mpi4py import MPI


def plot_history(z_array_reion_allmodels, 
                 lookback_array_reion_allmodels, cosmology_allmodels,
                 t_bigbang_allmodels, mass_frac_allmodels, 
                 model_tags, output_dir, output_tag, output_format,
                 passed_ax=None):
    """
    Plots the evolution of the neutral hydrogen fraction as a function of
    time since Big Bang.

    .. note::
        The snapshot range only covers reionization (not the full simulation
        snapshot range).

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

    mass_frac_allmodels : 2D nested list of floats. Dimensions equal to
                          ``z_array_reion_allmodels``.
        The mass weighted neutral fraction at each snapshot for each model.

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    passed_ax : ``matplotlib`` axes, optional
        If defined, the history will be plotted onto the passed axis and
        returned.  The figure will not be saved. 

    Returns
    ---------

    If ``passed_ax = None``, None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".

    If ``passed_ax`` is passed, the axis is returned with the ionization
    history plotted.
    """
    if passed_ax:
        ax1 = passed_ax
    else:
        fig1 = plt.figure(figsize = (8,8))
        ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_reion_allmodels)):

        ax1.plot(lookback_array_reion_allmodels[model_number], 
                 mass_frac_allmodels[model_number], 
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

    ax2 = ps.add_time_z_axis(ax1, cosmology_allmodels[0],
                             t_bigbang_allmodels[0]/1.0e3)

    ax1 = ps.adjust_axis(ax1, ps.global_axiswidth, ps.global_tickwidth,
                         ps.global_major_ticklength, ps.global_minor_ticklength) 

    if not passed_ax:
        leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(ps.global_legendsize)

    if passed_ax:
        return ax1
    else:
        outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
        fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile1))
        plt.close(fig1)

        return None


def plot_nion(z_array_reion_allmodels, lookback_array_reion_allmodels,
              cosmology_allmodels, t_bigbang_allmodels, nion_allmodels,
              model_tags, output_dir, output_tag, output_format):
    """
    Plots the ionizing emissivity as a function of time since Big Bang. 

    .. note::
        The snapshot range only covers reionization (not the full simulation
        snapshot range).

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
        The ionizing emissivity at each snapshot for each model. Units is
        1.0e50 Photons/Mpc^3. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_tag : String.
        Tag added to the name of the output file.

    output_dir : String
        Directory where the plot is saved.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_reion_allmodels)):

        # Units for photons are 1.0e50 so move to log properly.
        log_sum = 50 + np.log10(nion_allmodels[model_number])

        ax1.plot(lookback_array_reion_allmodels[model_number], 
                 log_sum, 
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

    ax1 = ps.plot_bouwens2015(cosmology_allmodels[0], t_bigbang_allmodels[0], ax1)

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

    ax1 = ps.adjust_axis(ax1, ps.global_axiswidth, ps.global_tickwidth,
                         ps.global_major_ticklength, ps.global_minor_ticklength) 

    leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


def plot_ps_fixed_XHI(k, P21, PHII, fixed_XHI_values, model_tags, output_dir,
                      output_tag, output_format):
    """
    Plots the 21cm power spectrum at specified fixed neutral fractions. 

    Parameters
    ----------

    k : 3D nested list of floats. Outer length is number of models. Next is
        number of HI fractions we're plotting (``fixed_XHI_values``) then 
        finally is number of bins.
        The wavenumber the spectra are binned on. Units are h/Mpc. 

    P21, PHII : 3D nested list of floats. Dimensions are identical to ``k``.
        The 21cm and HII power spectra. 21cm units are mK and HII is unitless.

    fixed_XHI_values : List of floats
        The mass-averaged neutral hydrogen fraction we're plotting at.

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig, ax = plt.subplots(nrows = 1, ncols = len(fixed_XHI_values),
                           sharey='row', figsize=(16, 6))

    for model_number in range(len(k)):
        for fraction in range(len(fixed_XHI_values)):

            this_ax = ax[fraction]

            w = np.where(k[model_number][fraction] > 9e-2)[0]
            label = model_tags[model_number]
 
            this_ax.plot(k[model_number][fraction][w],
                         P21[model_number][fraction][w],
                         color = ps.colors[model_number],
                         ls = ps.linestyles[model_number],
                         lw = 2, rasterized=True, label = label)
            
            if model_number != 0:
                continue

            this_ax = ps.adjust_axis(this_ax, ps.global_axiswidth,
                                     ps.global_tickwidth,
                                     ps.global_major_ticklength,
                                     ps.global_minor_ticklength)


            this_ax.set_xlabel(r'$\mathbf{k \: \left[Mpc^{-1}h\right]}$',
                               size = ps.global_labelsize)
            this_ax.set_xscale('log')
            this_ax.set_xlim([7e-2, 6])

            HI_string = "{0:.2f}".format(fixed_XHI_values[fraction])

            label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"
            this_ax.text(0.05, 0.9, label, transform = this_ax.transAxes,
                         fontsize = ps.global_fontsize)


            #tick_locs = np.arange(-3, 2.5, 1.0)
            #this_ax.set_xticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
            #                        fontsize = ps.global_fontsize)

    ax[0].set_ylabel(r'$\mathbf{\Delta_{21}^2 \left[mK^2\right]}$',
                     size = ps.global_labelsize)
    ax[0].set_yscale('log', nonposy='clip')

    #tick_locs = np.arange(-1.0, 5.5, 1.0)
    #ax[0].set_yticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
    #                      fontsize = ps.global_fontsize)

    leg = ax[0].legend(loc='lower right', numpoints=1,
                       labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize-2)

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)


    outputFile = "{0}/{1}.{2}".format(output_dir,
                                     output_tag,
                                     output_format)
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


def plot_incomp_contour(reion_comp_allmodels, alpha, beta, ax):
    """
    For some models of fesc reionization has not completed. When we plot
    contours (of constant tau or duration), we wish to mark which parameter
    choices did not complete reionization.

    This function plots a shaded regions to cover these incomplete models. 

    Parameters
    ----------

    reion_comp_allmodels : List of integers. Length is number of models. 
        Flag which denotes whether reionization has been completed for a given
        model.

    alpha : List of floats. 
        Range of alpha values used for the fesc models. 

    beta : List of floats. 
        Range of beta values used for the fesc models. 

    ax : ``matplotlib`` axis object.
        The axis we're plotting the shaded region on.

    Returns
    ---------

    ax : ``matplotlib`` axis object.
        The axis with the shaded region plotted on it.
    """

    incomp_alpha_inds = np.where(reion_comp_allmodels[0] == 0)[0]
    incomp_beta_inds = np.where(reion_comp_allmodels[1] == 0)[0]

    if alpha[0] == beta[0] == 0:
        incomp_alpha_low = 0.0
        incomp_beta_low = 0.0
    else:
        incomp_alpha_low = incomp_alpha_inds[0] * (alpha[1] - alpha[0]) 
        incomp_beta_low = incomp_beta_inds[0] * (beta[1] - beta[0]) 

    incomp_alpha_high = incomp_alpha_inds[-1] * (alpha[1] - alpha[0]) 
    incomp_beta_high = incomp_beta_inds[-1] * (beta[1] - beta[0]) 

    print("Smallest (alpha, beta) not completed is ({0}, "
          "{1})".format(incomp_alpha_low, incomp_beta_low))
    print("Largest (alpha, beta) not completed is ({0}, "
          "{1})".format(incomp_alpha_high, incomp_beta_high))

    line_low = -alpha
    line_high = -(0.1/0.3)*alpha + 0.1

    ax.fill_between(alpha, line_low, line_high, alpha = 0.3)

    return ax


def plot_tau_contours(tau_highz, reion_completed, alpha_beta_limits,
                      output_dir, output_tag, output_format):
    """
    Plots contours of constant optical depth tau. We recommend only plotting
    this if you're read the data for different alpha/beta normalization values. 

    Parameters
    ----------

    tau_highz: List of floats. Length is number of models. 
        The optical depth at the highest redshift (snap 0).

    reion_completed : List of integers. Length is number of models.
        Flag which denotes whether reionization has been completed for a given
        model.

    alpha_beta_limits : 2x3 list of floats.
        Defines the minimum, maximum and step size of the alpha and beta values
        that we read data from.  The 0th row is the alpha values, 1st is beta.

    output_tag : String.
        Tag added to the name of the output file.

    output_dir : String
        Directory where the plot is saved.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    tau_highz= np.array(tau_highz)
    reion_completed = np.array(reion_completed)

    # Need to add the (alpha, beta) = (0.0, 0.0) model because it has no
    # reionization (tau is not defined).
    if alpha_beta_limits[0][0] == alpha_beta_limits[1][0] == 0:
        tau_highz = np.insert(tau_highz, [0], np.nan)
        reion_completed = np.insert(reion_completed, [0], 0)

    alpha_low = alpha_beta_limits[0][0]
    alpha_high = alpha_beta_limits[0][1]
    alpha_step = alpha_beta_limits[0][2] 
    alpha = np.arange(alpha_low, alpha_high + alpha_step, alpha_step)

    beta_low = alpha_beta_limits[1][0]
    beta_high = alpha_beta_limits[1][1]
    beta_step = alpha_beta_limits[1][2]
    beta = np.arange(beta_low, beta_high + beta_step, beta_step)

    X, Y = np.meshgrid(beta, alpha)
    tau = tau_highz.reshape(len(X), len(Y))
    reion_completed = reion_completed.reshape(len(X), len(Y))

    CS = ax1.contour(Y, X, tau)

    ax1 = plot_incomp_contour(reion_completed, alpha, beta, ax1)

    ax1.clabel(CS, inline=1, fontsize = ps.global_fontsize)

    ax1.set_xlabel(r"$\mathbf{\alpha}$", size = ps.global_fontsize)
    ax1.set_ylabel(r"$\mathbf{\beta}$", size = ps.global_fontsize)

    ax1 = ps.adjust_axis(ax1, ps.global_axiswidth, ps.global_tickwidth,
                         ps.global_major_ticklength, ps.global_minor_ticklength) 

    ax1.set_xlim([alpha_low, alpha_high])
    ax1.set_ylim([beta_low, beta_high])

    ax1.xaxis.set_major_locator(mtick.MultipleLocator(alpha_step))
    ax1.yaxis.set_major_locator(mtick.MultipleLocator(beta_step))

    tick_locs = np.arange(alpha_low - alpha_step, alpha_high + alpha_step,
                          alpha_step) 
    ax1.set_xticklabels([r"$\mathbf{%.1f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)

    tick_locs = np.arange(beta_low - beta_step, beta_high + beta_step,
                          beta_step) 
    ax1.set_yticklabels([r"$\mathbf{%.2f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)

    outputFile = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


def plot_duration_contours(duration_t, reion_completed, alpha_beta_limits,
                           output_dir, output_tag, output_format):
    """
    Plots contours of constant reionization duration. We recommend only plotting
    this if you're read the data for different alpha/beta normalization values. 

    Parameters
    ----------

    duration_t : 2D nested list of floats. Outer length is number of models,
                 inner is 3.
        The lookback time corresponding to the start, middle and end of
        reionization. Units is Myr.

    reion_completed : List of integers. Length is number of models.
        Flag which denotes whether reionization has been completed for a given
        model.

    alpha_beta_limits : 2x3 list of floats.
        Defines the minimum, maximum and step size of the alpha and beta values
        that we read data from.  The 0th row is the alpha values, 1st is beta.

    output_tag : String.
        Tag added to the name of the output file.

    output_dir : String
        Directory where the plot is saved.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    dt = np.zeros(len(duration_t))
    for model_number in range(len(duration_t)):
        dt[model_number] = duration_t[model_number][-1] - \
                           duration_t[model_number][0]

    # First generate the duration of reionization
    if alpha_beta_limits[0][0] == alpha_beta_limits[1][0] == 0:
        reion_completed = np.insert(reion_completed, [0], 0)
        dt = np.insert(dt, [0], np.nan)

    alpha_low = alpha_beta_limits[0][0]
    alpha_high = alpha_beta_limits[0][1]
    alpha_step = alpha_beta_limits[0][2] 
    alpha = np.arange(alpha_low, alpha_high + alpha_step, alpha_step)

    beta_low = alpha_beta_limits[1][0]
    beta_high = alpha_beta_limits[1][1]
    beta_step = alpha_beta_limits[1][2]
    beta = np.arange(beta_low, beta_high + beta_step, beta_step)

    X, Y = np.meshgrid(beta, alpha)
    dt = dt.reshape(len(X), len(Y))
    reion_completed = reion_completed.reshape(len(X), len(Y))

    CS = ax1.contour(Y, X, dt)

    ax1 = plot_incomp_contour(reion_completed, alpha, beta, ax1)    

    ax1.clabel(CS, inline=1, fontsize = ps.global_fontsize)

    ax1.set_xlabel(r"$\mathbf{\alpha}$", size = ps.global_fontsize)
    ax1.set_ylabel(r"$\mathbf{\beta}$", size = ps.global_fontsize)

    ax1 = ps.adjust_axis(ax1, ps.global_axiswidth, ps.global_tickwidth,
                         ps.global_major_ticklength, ps.global_minor_ticklength) 

    ax1.set_xlim([alpha[0], alpha[-1]])
    ax1.set_ylim([beta[0], beta[-1]])

    ax1.xaxis.set_major_locator(mtick.MultipleLocator(0.1))
    ax1.yaxis.set_major_locator(mtick.MultipleLocator(0.05))

    tick_locs = np.arange(min(alpha)-0.1, max(alpha)+0.1, 0.1)
    ax1.set_xticklabels([r"$\mathbf{%.1f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)

    tick_locs = np.arange(min(beta)-0.05, max(beta)+0.05, 0.05)
    ax1.set_yticklabels([r"$\mathbf{%.2f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)

    outputFile = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


def plot_tau(z_array_reion_allmodels, lookback_array_reion_allmodels,
             cosmology_allmodels, t_bigbang_allmodels, tau_allmodels,
             model_tags, output_dir, output_tag, output_format,
             passed_ax=None):
    """
    Plots the evolution of the Thomson optical depth tau as a function of
    time since Big Bang.

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

    tau_allmodels : 2D nested list of floats. Dimensions equal to 
                    ``z_array_reion_allmodels``.
        The Thomson optical depth for each model at each snapshot.

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    passed_ax : ``matplotlib`` axes, optional
        If defined, the history will be plotted onto the passed axis and
        returned.  The figure will not be saved. 

    Returns
    ---------

    If ``passed_ax = None``, None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".

    If ``passed_ax`` is passed, the axis is returned with the ionization
    history plotted.
    """
   
 
    if passed_ax:
        ax1 = passed_ax
    else:
        fig1 = plt.figure(figsize = (8,8))
        ax1 = fig1.add_subplot(111)

    for model_number in range(len(tau_allmodels)):
        ax1.plot(lookback_array_reion_allmodels[model_number], 
                 tau_allmodels[model_number], 
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])
        
    ax1.set_xlabel(r"$\mathbf{Time \: since \: Big \: Bang \: [Myr}]$",
                  size = ps.global_labelsize)

    ax1.set_ylabel(r"$\mathbf{\tau}$",
                 size = ps.global_labelsize)
    ax1.set_ylim([0.042, 0.072])
    tick_locs = np.arange(0.04, 0.075, 0.005)
    ax1.set_yticklabels([r"$\mathbf{%.3f}$" % x for x in tick_locs],
                        fontsize = ps.global_fontsize)

    ax1 = ps.adjust_axis(ax1, ps.global_axiswidth, ps.global_tickwidth,
                         ps.global_major_ticklength, ps.global_minor_ticklength) 
    ax2 = ps.add_time_z_axis(ax1, cosmology_allmodels[0],
                             t_bigbang_allmodels[0]/1.0e3)

    ax1.fill_between(np.arange(200, 1200, 0.01), 0.058 - 0.012, 0.058 + 0.012,
                     color = 'k', alpha = 0.3)
    ax1.axhline(y = 0.058, xmin = 0, xmax = 20, color = 'k', alpha = 0.3)
    ax1.text(850, 0.0570, r"$\mathrm{Planck \: 2016}$")

    if not passed_ax:
        leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(ps.global_legendsize)
 
    if passed_ax:
        return ax1
    else:
        outputFile = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()


def plot_combined_history_tau(z_array_reion_allmodels,
                              lookback_array_reion_allmodels,                              
                              cosmology_allmodels, t_bigbang_allmodels, 
                              mass_frac_allmodels, tau_allmodels, model_tags,
                              output_dir, output_tag, output_format):
    """
    Plots the evolution of ionized hydrogen and Thomson optical depth (as a
    function of time since Big Bang) on a single plot.

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

    tau_allmodels : 2D nested list of floats. Dimensions equal to 
                    ``z_array_reion_allmodels``.
        The Thomson optical depth for each model at each snapshot.

    mass_frac_allmodels : 2D nested list of floats. Dimensions equal to
                          ``z_array_reion_allmodels``.
        The mass weighted neutral fraction at each snapshot for each model.

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    passed_ax : ``matplotlib`` axes, optional
        If defined, the history will be plotted onto the passed axis and
        returned.  The figure will not be saved. 

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    fig, ax = plt.subplots(nrows = 1, ncols = 2, 
                           sharex=False, sharey=False, figsize=(16, 8))

    # Pass each axis to the correct function to create the plots.
    ax[0] = plot_history(z_array_reion_allmodels, 
                         lookback_array_reion_allmodels, cosmology_allmodels,
                         t_bigbang_allmodels, mass_frac_allmodels,
                         model_tags, output_dir, output_tag, output_format,
                         passed_ax=ax[0])

    ax[1] = plot_tau(z_array_reion_allmodels, lookback_array_reion_allmodels,
                     cosmology_allmodels, t_bigbang_allmodels, tau_allmodels,
                     model_tags, output_dir, output_tag, output_format,
                     passed_ax=ax[1])

    leg = ax[0].legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    plt.tight_layout()    

    outputFile = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


def plot_ps_scales(k_allmodels, P21_allmodels, PHII_allmodels,
                   mass_frac_allmodels, fixed_XHI_values,
                   small_scale_def, large_scale_def, model_tags, output_dir,
                   output_tag, output_format):
    """
    Plots the small scale 21cm power spectrum as a function of the large scale
    21cm power. 

    We also mark the specified ``fixed_XHI_values`` to draw the eye.

    .. note::
        The snapshot range only covers reionization (not the full simulation
        snapshot range).

    Parameters
    ----------

    k_allmodels : 3D nested list of floats. Outer length is number of models.
                  Next is the number of snapshots then finally is the number of
                  bins. 
        The wavenumber the spectra are binned on. Units are h/Mpc.

    P21, PHII : 3D nested list of floats. Dimensions are identical to ``k``.
        The 21cm and HII power spectra. 21cm units are mK and HII is unitless.

    mass_frac_allmodels : 2D nested list of floats. Outer length is number of
                          models, inner is number of snapshots. 
        The mass weighted neutral fraction at each snapshot for each model.

    fixed_XHI_values : List of floats
        The mass-averaged neutral hydrogen fraction we're marking on the plot. 

    small_scale_def, large_scale_def : Floats.
        The wavenumber values (in h/Mpc) that correspond to 'small' and 'large'
        scales.

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    num_models = len(k_allmodels)

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    snap_idx_target = []

    k_small_scale = []
    k_large_scale = []

    P21_small_scale = []
    P21_large_scale = []

    PHII_small_scale = []
    PHII_large_scale = []

    # We first need to find the power on the specified scales.
    for model_number in range(num_models):
        snap_idx_target.append([])

        k_small_scale.append([])
        k_large_scale.append([])

        P21_small_scale.append([])
        P21_large_scale.append([])

        PHII_small_scale.append([])
        PHII_large_scale.append([]) 

        k_this_model = k_allmodels[model_number]
        P21_this_model = P21_allmodels[model_number]
        PHII_this_model = PHII_allmodels[model_number]

        # For all the snapshots find the values at the specified scales.
        for snap_idx in range(len(k_this_model)):
            small_idx = (np.abs(k_this_model[snap_idx] - small_scale_def)).argmin()
            large_idx = (np.abs(k_this_model[snap_idx] - large_scale_def)).argmin()
            
            k_small_scale[model_number].append(k_this_model[snap_idx][small_idx])
            P21_small_scale[model_number].append(P21_this_model[snap_idx][small_idx])
            PHII_small_scale[model_number].append(PHII_this_model[snap_idx][small_idx])

            k_large_scale[model_number].append(k_this_model[snap_idx][large_idx])
            P21_large_scale[model_number].append(P21_this_model[snap_idx][large_idx])
            PHII_large_scale[model_number].append(PHII_this_model[snap_idx][large_idx])

        # Now we have the power we can plot!
        ax1.plot(P21_small_scale[model_number],
                 P21_large_scale[model_number],
                 color = ps.colors[model_number],
                 ls = '-',
                 label = model_tags[model_number],
                 zorder = 1) 

        # Find the index corresponding to the special HI fractions.
        for count, val in enumerate(fixed_XHI_values):
            idx = (np.abs(mass_frac_allmodels[model_number] - val)).argmin()
            snap_idx_target[model_number].append(idx)

            ax1.scatter(P21_small_scale[model_number][idx],
                        P21_large_scale[model_number][idx],
                        s = 100, rasterized = True,
                        marker = ps.markers[count],
                        color = ps.colors[model_number], 
                        linewidths = 2,
                        zorder = 2)

    # Make a one-to-one line and plot it.
    one_to_one = np.arange(0.0, 40.0, 1e-6) 
    ax1.plot(one_to_one, one_to_one, ls = '--', lw = 0.5, color = 'k')

    small_scale_string = "{0:.1f}".format(small_scale_def)
    x_label = r"$\mathbf{\Delta_{21}^2 \left[mK^2\right] (k = " + \
              small_scale_string + r"Mpc^{-1}h)}$" 
    ax1.set_xlabel(x_label, size = ps.global_labelsize)

    large_scale_string = "{0:.1f}".format(large_scale_def)
    y_label = r"$\mathbf{\Delta_{21}^2 \left[mK^2\right] (k = " + \
              large_scale_string + r"Mpc^{-1}h)}$" 
    ax1.set_ylabel(y_label, size = ps.global_labelsize)

    max_smallscale = np.ceil(np.max(P21_small_scale))
    max_small_quot, max_small_rem = divmod(max_smallscale, 5) 

    ax1.set_xlim([0.0, (max_small_quot+1)*5])
    ax1.set_ylim([0.0, (max_small_quot+1)*5]) 

    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(1))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(1))

    ax1.xaxis.set_major_locator(mtick.MultipleLocator(5))
    ax1.yaxis.set_major_locator(mtick.MultipleLocator(5))

    ax1 = ps.adjust_axis(ax1, ps.global_axiswidth, ps.global_tickwidth,
                         ps.global_major_ticklength, ps.global_minor_ticklength) 

    for count, val in enumerate(fixed_XHI_values):
        HI_string = "{0:.2f}".format(val)
        label = r"$\mathbf{\langle \chi_{HI}\rangle = " + HI_string + r"}$"
        ax1.scatter(np.nan, np.nan,
                    s = 60, rasterized = True,
                    marker = ps.markers[count],
                    color = "None",
                    edgecolors = 'k', 
                    linewidths = 2,
                    zorder = 2, label = label)

    leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1,
                     markerscale=6)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

def plot_slices_XHI(z_array_reion_allmodels, cosmology_allmodels,
                    mass_frac_allmodels, XHII_fbase_allmodels,
                    XHII_precision_allmodels, GridSize_allmodels,
                    boxsize_allmodels, first_snap_allmodels,
                    fixed_XHI_values, cut_slice, cut_thickness, model_tags,
                    output_dir, output_tag, output_format):
    """
    Plots slices of the ionization fields for each model at fixed neutral
    hydrogen fractions.

    .. note::
        The snapshot range only covers reionization (not the full simulation
        snapshot range).

    Parameters
    ----------

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    cosmology_allmodels : List of class ``astropy.cosmology``. Length is number
                          of models.
        ``Astropy`` class containing the cosmology for each model.

    mass_frac_allmodels : 2D nested list of floats. Dimensions equal to
                          ``z_array_reion_allmodels``.
        The mass weighted neutral fraction at each snapshot for each model.

    XHII_fbase_allmodels : List of strings. Length is number of models.
        The base filename for the ionization fields for each model.

    XHII_precision_allmodels : List of integers. Length is number of models.
        The precision of the ionization fields for each model.
        1 : Float, 2 : Double.

    GridSize_allmodels : List of integers. Length is number of models.
        The number of grid cells (along a box size) for each model.

    boxsize_allmodels : List of integers. Length is number of models.
        The simulation box size for each model (units are Mpc/h).

    first_snap_allmodels : List of integers. Length is number of models.
        The snapshot where ``cifog`` starts calculations for each model.

    fixed_XHI_values : List of floats.
        The neutral hydrogen fractions we're calculating the slices at.
        Defined by the user in ``paper_plots.py``.

    cut_slice : Integer.
        The grid index for model 0 that we're cutting at. All other models are
        normalized to ensure the cut is at the same spatial position at model 0.

    cut_thickness : Integer.
        The thickness (in grid cells) that the model 0 slice will be. All other
        models normalized to ensure the same spatial thickness as model 0. 

    model_tags : List of strings. Length is number of models.
        Legend entry for each model.

    output_dir : String
        Directory where the plot is saved.

    output_tag : String.
        Tag added to the name of the output file.

    output_format : String
        Format the plot is saved in.

    Returns
    ---------

    None. The figure is saved as "<output_dir>/<output_tag>.<output_format>".
    """

    num_models = len(mass_frac_allmodels)
    num_fractions = len(fixed_XHI_values)

    fig, ax = plt.subplots(nrows=num_models, ncols=num_fractions,
                           sharey=False, sharex=False, figsize=(12, 12))
    fig.subplots_adjust(wspace = 0.02, hspace = 0.02)

    # Each model can have a different gridsize.  We want to cut at the same
    # spatial location for each model so we will normalize to model 0.
    mod0_gridsize = GridSize_allmodels[0]
    mod0_boxsize = boxsize_allmodels[0]

    for model_number in range(num_models):
        model_gridsize = GridSize_allmodels[model_number]
        model_boxsize = boxsize_allmodels[model_number]

        for frac_number, frac_val in enumerate(fixed_XHI_values):

            this_ax = ax[frac_number, model_number]

            # First find the snapshot that corresponds to this XHI value.
            snap_idx = (np.abs(mass_frac_allmodels[model_number] - frac_val)).argmin()
            snap_z = z_array_reion_allmodels[model_number][snap_idx]

            # These are cifog files so need to add 1.
            cifog_snapnum = snap_idx + first_snap_allmodels[model_number] + 1

            XHII_path = "{0}_{1:03d}".format(XHII_fbase_allmodels[model_number],
                                             cifog_snapnum)
            XHII = rs.read_binary_grid(XHII_path, GridSize_allmodels[model_number],
                                       XHII_precision_allmodels[model_number])

            # Set this up to get nice plotting.
            ionized_cells = np.log10(1 - XHII)

            # Find the grid index that corresponds to the same spatial scale as
            # model 0.
            index_cut = int(cut_slice *  model_boxsize/mod0_boxsize * \
                            model_gridsize / mod0_gridsize)
            # Then we use the `thickness cut` variable of model 0 and then
            # scale other models so the same thickness is being cut.
            thickness_cut = int(np.ceil(cut_thickness * model_boxsize/mod0_boxsize * \
                                        model_gridsize / mod0_gridsize))

            im = this_ax.imshow(ionized_cells[:,:,index_cut:index_cut+thickness_cut].mean(axis=-1),
                                interpolation="none", origin="low",
                                extent = [0.0, model_boxsize, 0.0, model_boxsize],
                                vmin=-8, vmax=0, cmap="afmhot_r")
            this_ax.axis('off')

            this_ax.set_xlim([0.0, model_boxsize])
            this_ax.set_ylim([0.0, model_boxsize])

            # The fraction strings are plotted on the left-side. So set them
            # for model 0 (column 0).
            if model_number == 0:
                HI_string = "{0:.2f}".format(frac_val)

                label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"
                this_ax.text(-0.2,0.8, label, transform=this_ax.transAxes, 
                             size=ps.global_labelsize - 10,
                             rotation=90)
            # The model strings are plotted on the top. So set them for
            # fraction 0 (row 0).
            if frac_number == 0:
                title = model_tags[model_number]                 
                this_ax.set_title(title, 
                                  size = ps.global_labelsize - 4) 

            # Finally add the redshift to the top right of the axis.
            z_label = r"$z = %.2f$" %(snap_z)
            z_text = this_ax.text(0.55,0.9, z_label,
                                  transform=this_ax.transAxes, 
                                  size=ps.global_labelsize - 16,
                                  color='k')
            # Add a white background to make the label legible.
            z_text.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
            plt.draw()

        # HI Fraction Loop.
    # Model Loop.

    # All the models have been plotted. Now lets fix up the colorbar.
    cax = fig.add_axes([0.91, 0.11, 0.03, 0.77])
    ticks = np.arange(-8.0, 1.0, 1.0)
    cbar = fig.colorbar(im, cax=cax, ticks=ticks) 
    cbar.ax.set_yticklabels([r"$\mathbf{%d}$" % x for x in ticks], 
                            fontsize = ps.global_legendsize+10)
    cbar.ax.set_ylabel(r'$\mathbf{log_{10}\left(\chi_{HI}\right)}$',
                       rotation = 90, size = ps.global_labelsize)
    #cbar.ax.tick_params(labelsize = ps.global_legendsize + 10)

    # Done! Save time.
    outputFile = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

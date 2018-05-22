#!/usr/bin/env python 
import matplotlib
matplotlib.use('Agg') 
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

from astropy import units as u
from astropy import cosmology
from matplotlib import gridspec

import PlotScripts
import ReadScripts
import AllVars

m_low = 7.0 
m_high = 15.0

m_gal_low = 5.0
m_gal_high = 15.0

bin_width = 0.2

def my_load_data(fname):

    data = ReadScripts.load_data(fname)

    fesc = data[:,0]
    fej = data[:,1]  # Fraction of baryons in the ejected reservoir.
    Len = data[:,2]  # Number of particles in the host halo. 
    halo_mass= data[:,3]  # Virial mass of host halo (1.0e10 Msun/h). 
    fof_halo_mass = data[:,4]  # Virial mass of background FoF Halo
                               # (1.0e10Msun/h). 
    stellar_mass= data[:,5]  # Stellar mass of galaxy (1.0e10 Msun/h). 
    Ngamma_HI = data[:,6]  # Number of emitted HI ionizing photons
                           # 1.0e50 photon/s. 

    halo_mass = halo_mass * 1.0e10 / AllVars.Hubble_h
    fof_halo_mass = fof_halo_mass * 1.0e10 / AllVars.Hubble_h
    stellar_mass = stellar_mass * 1.0e10 / AllVars.Hubble_h
    Ngamma_HI = Ngamma_HI * 1.0e50

    return fesc, fej, Len, halo_mass, fof_halo_mass, stellar_mass, Ngamma_HI


def plot_fesc_z(filebase, PlotSnapshot, model_tags, output_tag,
                output_format=".png", my_halo_mlow=m_low, my_halo_mhigh=m_high,
                my_gal_mlow=m_gal_low, my_gal_mhigh=m_gal_high):

    def adjust_halo_plot(ax):
        # Place horizontal lines for the quasar prescription values.
        ax.axhline(0.20, 0, 100, color ='k', 
                    lw = PlotScripts.global_linewidth, ls = '--')
        ax.axhline(0.10, 0, 100, color ='k', 
                    lw = PlotScripts.global_linewidth, ls= '-.')
        ax.text(11.0, 0.22, r"$f_\mathrm{esc} = 0.20$", color = 'k', 
                size = PlotScripts.global_fontsize)
        ax.text(9.25, 0.11, r"$f_\mathrm{esc, \: base}$", color = 'k', 
                size = PlotScripts.global_fontsize)
 
        ax.set_xlabel(r'$\log_{10}\: \mathrm{M}_\mathrm{vir}\ [\mathrm{M}_{\odot}]$', 
                      size = PlotScripts.global_labelsize)
        ax.set_ylabel(r'$\langle f_\mathrm{esc}\rangle_\mathrm{M_{vir}}$', 
                      size = PlotScripts.global_labelsize)
        leg = ax.legend(loc='upper left', bbox_to_anchor=(0.3, 1.02), 
                        numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)

    fig2 = plt.figure(figsize=(8,8))
    ax3 = fig2.add_subplot(111)

    for model_number in range(len(filebase)):
        for count, snap in enumerate(PlotSnapshot[model_number]):
            fname = "{0}_{1:03d}".format(filebase[model_number], snap)

            fesc, fej, Len, halo_mass, fof_halo_mass,\
            stellar_mass, Ngamma_HI = my_load_data(fname)

            w = np.where(stellar_mass > 0.0)[0]

            Mvir = np.log10(halo_mass[w])
            Mstar = np.log10(stellar_mass[w])
            print("There are {0} galaxies at Snapshot {1}".format(len(w), 
                                                                  snap))

            mean_fescMvir, std_fescMvir, N_fescMvir, _, bins_mid_Mvir = \
            AllVars.Calculate_2D_Mean(Mvir, fesc[w], bin_width, 
                                      my_halo_mlow, my_halo_mhigh)

            mean_fescMstar, std_fescMstar, N_fescMstar, _, bins_mid_Mstar = \
            AllVars.Calculate_2D_Mean(Mstar, fesc[w], bin_width, 
                                      my_gal_mlow, my_gal_mhigh)


            if model_number == 0:
                label = r"$\mathbf{z = " + str(int(round(AllVars.SnapZ[snap]))) + "}$"
            else:
                label = ""

            w_halo = np.where(N_fescMvir > 0.0)[0]
            ax1.plot(bins_mid_Mvir[w_halo], mean_fescMvir[w_halo],
                     color = PlotScripts.colors[count], 
                     ls = PlotScripts.linestyles[model_number],
                     label = label, lw = PlotScripts.global_linewidth)

    adjust_halo_plot(ax1)

    outputFile1 = './{0}_halo{1}'.format(output_tag, output_format)
    outputFile2 = './{0}_galaxy{1}'.format(output_tag, output_format)

    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure

    print('Saved file to {0}'.format(outputFile1))
    print('Saved file to {0}'.format(outputFile2))
    plt.close()


if __name__ == '__main__':

    AllVars.Set_Params_Kali()
    PlotScripts.Set_Params_Plot()

    fname=["/lustre/projects/p004_swin/jseiler/kali/tracking/grids/nion/properties/test_misc_properties"]
    PlotSnapshot = [[40]]
    model_tags = ["test"]

    plot_fesc_z(fname, PlotSnapshot, model_tags, "test")


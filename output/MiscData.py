#!/usr/bin/env python
"""
Author: Jacob Seiler
Version: 0.1
"""

from __future__ import print_function

import os
from natsort import natsorted

def get_sage_ini_from_dir(directory, alpha_vals = None, beta_vals = None):

    SAGE_ini = []
    cifog_ini = []

    alpha_vals_strings = []
    if alpha_vals.any():
        for val in alpha_vals:
            string = "alpha{0:.1f}".format(val)
            alpha_vals_strings.append(string)

    beta_vals_strings = []
    if beta_vals.any():
        for val in beta_vals:
            string = "beta{0:.1f}".format(val)
            beta_vals_strings.append(string)

    for my_file in os.listdir(directory):
        if "SAGE" in my_file:
            if alpha_vals.any():
                if not any(val in my_file for val in alpha_vals_strings):
                    continue

            if beta_vals.any():
                if not any(val in my_file for val in beta_vals_strings):
                    continue
 
            fname = "{0}{1}".format(directory, my_file)
            SAGE_ini.append(fname)

        elif "cifog" in my_file:
            if alpha_vals.any():
                if not any(val in my_file for val in alpha_vals_strings):
                    continue

            if beta_vals.any():
                if not any(val in my_file for val in beta_vals_strings):
                    continue
 
            fname = "{0}{1}".format(directory, my_file)
            cifog_ini.append(fname)

    # Finally we want to sort the ini files in ascending order (using the
    # numbers within the string).
    SAGE_ini = natsorted(SAGE_ini)
    cifog_ini = natsorted(cifog_ini) 

    return SAGE_ini, cifog_ini

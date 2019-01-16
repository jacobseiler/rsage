"""
This function updates template ``SAGE`` and ``cifog`` ``.ini`` files to use the
specified directories. It should be used when using a brand new simulation as
it will allow the easy renaming and adjustment of the directories.

Note: This will NOT adjust any recipe flags or constants. It will only update
directory paths.

Author: Jacob Seiler
Version: 0.2.
"""

#!/usr/bin/env python
import os
import sys

# Get the directory the script is in.
# Used a global variable for convenience as quite a few functions will use this.
script_dir = os.path.dirname(os.path.realpath(__file__))
output_path = "{0}/output/".format(script_dir)
sys.path.append(output_path)

import ReadScripts as rs

utils_path = "{0}/utils".format(script_dir)
sys.path.append(utils_path)

import change_params as cp


def update_SAGE_dict(base_SAGE_params):
    """
    Wrapper function to update the parameters of the ``SAGE`` ``.ini`` file.

    Parameters
    ----------

    base_SAGE_params : Dictionary
        Dictionary keyed by the ``SAGE`` parameter field names and containing 
        the values from the ``.ini`` file. Refer to ``read_SAGE_ini`` in
        ``output/ReadScripts.py`` for the full key list.

    Returns
    ----------

    updated_dict : Dictionary
        Dictionary with identical data-structure to ``base_SAGE_params`` but
        with certain fields updated via user input.
    """

    # The fields that we want the user to update.
    fields_to_check = ["SimulationDir", "FileWithSnapList", "TreeName",
                       "TreeExtension", "FileNameGalaxies"]

    # Text to prompt the user input. 
    print_text = ["Path to the halo trees", 
                  "File with the simulation snapshot list",
                  "Prefix for the trees", "Suffix for the trees",
                  "Galaxy prefix name (used as unique run tag)"]

    # Some of these are required and others are optional.
    field_required = [False, False, False, False, True]
    check_type = ["dir", "file", "", "", ""]

    updated_dict = update_fields(base_SAGE_params, fields_to_check,
                                 print_text, field_required, check_type)

    return updated_dict


def update_cifog_dict(base_cifog_params):
    """
    Wrapper function to update the parameters of the ``cifog`` ``.ini`` file.

    Parameters
    ----------

    base_cifog_params : Dictionary
        Dictionary keyed by the ``cifog`` parameter field names and containing 
        the values from the ``.ini`` file. Refer to ``read_cifog_ini`` in
        ``output/ReadScripts.py`` for the full key list.

    Returns
    ----------

    updated_dict : Dictionary
        Dictionary with identical data-structure to ``base_cifog_params`` but
        with certain fields updated via user input.
    """

    # The fields that we want the user to update.
    fields_to_check = ["redshiftFile", "inputIgmDensityFile"]

    # Text to prompt user input.
    print_text = ["cifog redshift file", 
                  "Path to the input density file"]

    # How to we handle whether the update is required/optional.
    field_required = [False, False]
    check_type = ["file", ""]

    updated_dict = update_fields(base_cifog_params, fields_to_check,
                                 print_text, field_required, check_type)

    return updated_dict


def update_fields(base_dict, fields_to_check, print_text, field_required,
                  check_type):
    """
    For a subset of fields in a dictionary, prompts the user to provide updated
    values.

    Parameters
    ----------

    base_dict : Dictionary
        The fiducial dictionary we're updating the fields for. Will contain the
        ``.ini`` file parameters for either ``SAGE`` or ``cifog``.

    fields_to_check : List of strings
        The fields that we will be prompting the user to update.

    print_text : List of strings
        The text associated with the prompt for each field.

    field_required : List of boolean
        Specifies whether each of the fields in ``fields_to_check`` requires a
        user input. If ``False``, then the default value will be used.

    check_type : List of strings
        Some fields are paths to a directory or a file name. ``check_type``
        specifies which of these to check. Can be either ``"dir"`` or
        ``"file"``

    Returns
    ----------

    updated_dict : Dictionary
        Dictionary with identical data-structure to ``base_cifog_params`` but
        with certain fields updated via user input.
    """

    updated_dict = {}

    for idx in range(len(fields_to_check)):

        # Grab all the relevant info.
        field = fields_to_check[idx]
        text = print_text[idx]
        required = field_required[idx]
        check = check_type[idx] 

        if required:

            # If the field is marked as required, keep asking for a field until
            # the user enters one.
            my_field = None

            text_to_print = "{0} (must be given): ".format(text,
                                                           base_dict[field])
            while not my_field:
                my_field = input("{0}".format(text_to_print))
                if not my_field:
                    print("Must be specified.")

        else:
            # Otherwise just accept the default if they don't enter one.
            text_to_print = "{0} [default: {1}]: ".format(text,
                                                          base_dict[field])
            my_field = input("{0}".format(text_to_print))
            if not my_field: 
                my_field = base_dict[field]

        # Check that the directory path or the file actually exists.
        if check == "dir":
            if not os.path.isdir(my_field):
                print("Path {0} does not exist.".format(my_field))
                raise ValueError
        elif check == "file":
            if not os.path.isfile(my_field):
                print("File {0} does not exist.".format(my_field))
                raise ValueError

        updated_dict[field] = my_field

    return updated_dict


def adjust_ini():

    base_SAGE_ini = "{0}/ini_files/kali_SAGE.ini".format(script_dir)
    base_cifog_ini = "{0}/ini_files/kali_cifog.ini".format(script_dir)

    # Prompt for a SAGE ini path. If none provided, use the base.
    my_SAGE_ini = input("Template SAGE ini file [default: "
                        "{0}]: ".format(base_SAGE_ini))
    if not my_SAGE_ini:
        my_SAGE_ini = base_SAGE_ini

    # Do the same for cifog.
    my_cifog_ini = input("Template cifog ini file [default: "
                        "{0}]: ".format(base_cifog_ini))
    if not my_cifog_ini:
        my_cifog_ini = base_cifog_ini

    SAGE_params = rs.read_SAGE_ini(my_SAGE_ini)
    cifog_params, cifog_headers = rs.read_cifog_ini(my_cifog_ini)

    run_directory = None
    while not run_directory:
        run_directory = input("Base output directory: ")
        if not run_directory:
            print("Must be specified")

    SAGE_fields_update = update_SAGE_dict(SAGE_params)
    cifog_fields_update = update_cifog_dict(cifog_params)

    cp.create_directories(run_directory)
    cp.update_ini_files(base_SAGE_ini, base_cifog_ini,
                        SAGE_fields_update, cifog_fields_update,
                        run_directory)


if __name__ == '__main__':

    print("Welcome to the RSAGE ini file adjuster.")
    print("This script will adjust the directory names for both the SAGE and "
          "cifog ini files to be consistent with each other.") 
    print("For any of the following, if nothing is entered, we will use the "
          "default values for the template ini files specified.")
    print("")
    print("")

    adjust_ini()

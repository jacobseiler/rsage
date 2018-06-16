[![Build Status](https://travis-ci.org/jacobseiler/self_consistent_SAGE.svg?branch=master)](https://travis-ci.org/jacobseiler/self_consistent_SAGE)

# Description

This repository contains the Reionization using Semi-Analytic Galaxy Evolution (`RSAGE`) model.  This model is an augmented version of the Semi-Analytic Galaxy Evolution (`SAGE`) model (https://github.com/darrencroton/sage) that self-consistently couples galaxy evolution with the evolution of ionized gas during the Epoch of Reionization.  This coupling is achieved using the semi-numerical code [cifog](https://github.com/annehutter/grid-model). 

# Installation

## Pre-Requisites

`SAGE` requires only [GSL](https://www.gnu.org/software/gsl/) and should compile mostly out of the box.
`cifog` has slightly more complex build requirements and we refer to the `cifog` [README](https://github.com/annehutter/grid-model#pre-requisities) for full details on the required packages.

## Downloading 

```
$ git clone https://github.com/jacobseiler/rsage 
$ git submodule init
$ git submodule update 
``` 

`cifog` is kept as an independant repository and managed through the use of Git Submodules. To clone the version of `cifog` that is used within `RSAGE` (kept in its own separate [fork](https://github.com/jacobseiler/grid-model)), the `git submodule init` and `git submodule update` commands are used following the initial clone. 

## Building

`RSAGE` runs using three separate executables located in the *sage*, *grid-model* and *filter_mass*. 

```
$ cd sage/
$ make tests
$ cd grid-model/
$ make
$ cd filter_mass/
$ make 
```

# Running the code 

Like its parent model, `RSAGE` requires halo merger trees to run.  In addition to these trees, the computation of the ionization fields require dark matter overdensities grids. For testing purposes, we have included a small set of trees and dark matter density fields [here](BROKEN). 

Example parameter files for running the galaxy model and the semi-numerical reionization model are included in the *ini_files* directory (*kali.ini* and *cifog.ini* respectively).  For a complete description of the `cifog` parameters see the [README](https://github.com/annehutter/grid-model#parameter-file).  



# Authors and Reference Papers

The underlying semi-analytic galaxy evolution model was developed by [Darren Croton](https://github.com/darrencroton/sage) and is described in [Croton et al., 2016](https://arxiv.org/abs/1601.04709).
The semi-numerical model used to compute ionization fields was developed by [Anne Hutter](https://github.com/annehutter/grid-model) and is described in [Hutter, 2018](https://arxiv.org/abs/1803.00088).
The self-consistent coupling of these models was develpoed by Jacob Seiler and is described in (Stay tuned!). 

Any issues, questions or ranting can be sent to Jacob Seiler: jseiler@swin.edu.au 





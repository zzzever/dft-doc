---
title: Input Files
author: PWmat group
date: 2023-10-24
category: dft-doc
layout: post
---

PWmat needs a few basic input files to start the calculation: \hyperref[inputfile:etotinput]{parameter file} (must be named as {\bf etot.input}); \hyperref[inputfile:atomconfig]{structure file} (usually is {\bf atom.config} in our examples and tutorials); \hyperref[inputfile:pseudopotential]{pseudopotential files}. 

In some cases, one might also need to provide some \hyperref[optionalinput]{optional input files}, such as charge density (IN.RHO), high-symmetry-kpoints (IN.KPT), detailed solvent parameters (IN.SOLVENT). Some of them are simple so they can be written by hand. Some of them are binary files, which will be generated from the
previous calculations, then one can copy them to input file for the next calculation. For example, one should copy {\bf OUT.VR} to {\bf IN.VR} for non-self-consistent calculation. 

In the following, we will explain the long version of these files respectively.

## Parameter File

The parameter file is the most important input file for PWmat. It is a text file, and must be named as {\bf etot.input}. It contains all the parameters for the calculation. The parameters are divided into several groups, and each group is started with a keyword in the first line. The keywords are case-insensitive. The following is an example of parameter file:

```
JOB SCF
```
### JOB


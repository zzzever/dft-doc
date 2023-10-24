---
title: Input Files
author: PWmat group
date: 2023-10-24
category: dft-doc
layout: post
---

PWmat needs a few basic input files to start the calculation: `parameter file` (must be named as **etot.input**); `structure file` (usually is **atom.config** in our examples and tutorials); `pseudopotential files`. 

In some cases, one might also need to provide some `optional input files`, such as charge density (**IN.RHO**), high-symmetry-kpoints (**IN.KPT**), detailed solvent parameters (**IN.SOLVENT**). Some of them are simple so they can be written by hand. Some of them are binary files, which will be generated from the
previous calculations, then one can copy them to input file for the next calculation. For example, one should copy **OUT.VR** to **IN.VR** for non-self-consistent calculation. 

In the following, we will explain the long version of these files respectively.

# Parameter file (etot.input)
The parameter file must be named as **etot.input**. It is the most important input file, used to control how PWmat runs. Here is an example of the simplest:
```shell
4 1 
IN.ATOM = atom.config
IN.PSP1 = Si.NCPP.UPF
JOB = SCF
```
>TIP
>
>The first line must be two positive integers, which correspond to the tags **NODE1**, **NODE2** respectively. 
{: .block-tip }
    
The following lines in etot.input specify the name of the structure file, the name of the pseudopotential file and type of calculation.
You need to specify at least these parameters because they have no default values.

Except for the first line, the content has a format of **TAG** = **VALUE**. The orders of different tags can be arbitrarily changed. The names of the tags are case insensitive. One can add annotations after **\#** in line.
    
After running PWmat, one can also check the header of the output file **REPORT** and copy them as etot.input. 

In the following, we will explain the meaning of each tag.

## Control tags

### NODE1

### NODE2

### JOB

### ACCURACY

### PRECISION

### CONVERGENCE

### NUM_MPI_PER_GPU

### NUM_BLOCKED_PSI

### WF_STORE2DISK

### USE_GAUSSIAN

## System tags

### ECUT

### ECUT2

### ECUT2L

### ECUTP

### N123
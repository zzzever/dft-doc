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

```bash
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

| Tag | NODE1 |
| --- | --- |
| **Format** | the first integer in the first line |
| **Default** | none |
| **Related tags** | NODE2 |


The number of processors used to divide the G-space sphere and **N1 * N2 * N3** FFT grid. **N1 * N2 * N3** must be divisible by **NODE1**. 
The product of **NODE1** and **NODE2** is equal to the processors used by PWmat. 

>TIP
>
>If you run the same task with different **NODE1**, the automatically generated **N123** might be different, then the results might be slightly different.
{: .block-tip }

### NODE2

| Tag | NODE2 |
| --- | --- |
| **Format** | the second integer in the first line |
| **Default** | none |
| **Related tags** | NODE1 |

The number of processor groups to divide the k-points. One might check whether the number of k-points is divisible by **NODE2** for high efficiency. The larger **NODE2**, the larger the required memory.

The product of **NODE1** and **NODE2** is equal to the processors used by PWmat.

>WARNING
>
>**WARNING**: Hybrid functional calculation, k-point interpolation and electron-phonon coupling calculation do not support k-points parallelization. One must set **NODE2** = 1 in these cases.
{: .block-warning }

### JOB

| Tag | JOB |
| --- | --- |
| **Format** | JOB = [string] |
| **Default** | none |
| **Related tags** | none |

Controls what PWmat will do. **JOB** can be **SCF**, **NONSCF**, **DOS**, **MOMENT**, **RELAX**, **EGGFIT**, **MD**, **TDDFT**, **NAMD**, **NEB**, **DIMER**, **SCFEP**, **POTENTIAL**, **WKM**, **HPSI**, **ATOMIC\_ORB**, **TRANS**.

#### JOB = SCF
Do self-consistent field iterations.

|||
| --- | --- |
|**Related tags**| E\_ERROR, RHO\_ERROR, WG\_ERROR, FERMIDE, SCF\_ITER0\_*, CHARGE\_DECOMP, ... |
        
SCF calculates the charge density, the total energy. During SCF calculation, the atoms will not be moved.

Frequeuntly used `etot.input` settings for SCF calculation:

```bash
    4 1 
    IN.ATOM = atom.config 
    JOB = SCF 
    IN.PSP1 = Si.SG15.PBE.UPF 
    XCFUNCTIONAL = PBE 
    Ecut = 50 
    MP\_N123 = 9 9 9 0 0 0
```


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
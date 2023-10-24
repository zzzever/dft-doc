---
layout: home
title: PWmat Document
permalink: /
---

Make material simulation easy.

## Demo

Live demo on Github Pages: [https://zzzever.github.io](https://zzzever.github.io)

[![Static Badge](https://img.shields.io/badge/PWmat-doc-blue)](http://www.pwmat.com)



## Introduction

PWmat is a plane wave pseudopotential package for density functional
theory (DFT) calculations. It is designed to run efficiently on CPU/GPU processors. The best explanation of the algorithms used in PWmat
can be found in  [papers][1]. 

PWmat can perform the following calculations (set by the `JOB`
parameter in input file `etot.input`): 

|JOB|DESCRIPTION|
|---|---|
|SCF|self-consistent-field calculation|
|NONSCF|non-self-consistent-field calculations, e.g. for bandstructure calculation, which is usually done after a SCF calculation|
|DOS|density of state calculation, which is usually done after a NONSCF or SCF calculation, it is used to do partial density of state, or k-point interpolation|
|RELAX|atomic relaxation calculation and cell relaxation|
|EGGFIT|a preprocess fitting procedure to remove the egghead problem in RELAX|
|MD|ab initio molecular dynamics calculation|
|TDDFT|real-time time dependent density functional theory calculations|
|NAMD|non-adiabatic molecular dynamics, which is used to study carrier dynamics following a Born-Oppenheimer molecular dynamics|
|NEB|nudged elastic band calculation for barrier heights|
|DIMER|dimer method calculation for finding saddle points|
|POTENTIAL|generate DFT potential from input charge density|
|SCFEP|electron-phonon coupling constant calculation for a given pair of input electron states for all the phonon modes|
|WKM|a special Wannier Koopmann's method calculation for DFT band gap correction|
|ATOMIC\_ORB|atom's atomic wavefunction calculation|
|TRANS|quantum transport device calculation|

In order to run PWmat, one needs to provide the following necessary input files in the running
directory: 
1. `parameter file` (must be named as ***etot.input***); 
2. `structure file` (usually is ***atom.config***
in our examples and tutorials); 
3. `pseudopotential files`.

An example of workflow from the beginning to the end:
1. Get crystal structure from online database, or build it from visualization packages (e.g.[Q-Studio](https://mcloud.lonxun.com/){Q-Studio}).
2.  Download pseudopotential files from [PWmat website](http://www.pwmat.com/potential-download), then generate **etot.input** and convert structure file to PWmat format by using [PWkit](http://www.pwmat.com/pwkit-download).
3. Run command `mpirun -np 4 PWmat` to excute PWmat directly on a single computational node (Mstation). Or run `qsub job.pbs` if one is using TORQUE PBS on HPC clusters, `sbatch job.pbs` if one is using SLURM, to submit job on HPC clusters (like [Mcloud](https://mcloud.lonxun.com/)).
4. Collect the results by using [utilities](http://www.pwmat.com/utility-download). Prepare for the next calculation if needed.

All the above programs are pre-installed on the server. In the current PWmat releasing package, we also include one directory: `examples`,
which contains example cases for carrying out different jobs.

Besides to run PWmat alone, one may check our 
[modules](http://www.pwmat.com/module-download). 
They are recipes to carry out actual calculations for different scientific tasks. 
For example, calculating the defect level energies, the catalytic process in electrochemistry.
These are designed as tutorials to help users to finish the actual tasks. It can also be packages for more
sophisticated tasks, e.g., PyPWmat for phonon spectrum calculation, or using
YAMBO for GW calculation. We will develop more modules in the future, we also welcome our users to provide their own modules.

>##### TIP
> Due to rapid development, there could be some minor changes for
the tutorial, but the basic steps and ideas are the same, and the
changes should be obvious.

[1]: www.pwmat.com
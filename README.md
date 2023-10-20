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


[1]: www.pwmat.com
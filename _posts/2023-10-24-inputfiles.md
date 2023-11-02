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
    MP_N123 = 9 9 9 0 0 0
```

#### JOB = NONSCF
Do non-self-consistent calculations.

|||
| --- | --- |
|**Related tags**| IN.KPT, IN.VR, IN.NONSCF, ... |

NONSCF calculation requires an input potential, usually from a previous SCF calculation. One must set **IN.VR = T** in **etot.input**.

NONSCF calculates the eigen wave functions non-self-consistently, but do not calculate total energy. A denser K-mesh can be set, but for the bandstructure calculation, one can generate an explicitly K-path file and convert it by using **split\_kp.x**.

Frequeuntly used `etot.input` settings for NONSCF calculation:

```bash
    4 1 
    IN.ATOM = atom.config 
    JOB = NONSCF 
    IN.PSP1 = Si.SG15.PBE.UPF 
    XCFUNCTIONAL = PBE 
    Ecut = 50 
    IN.KPT = T   
    IN.VR = T   
```
In addition, some specific parameters can be set. Please refer to section [in.nonscf]() for details.

>
>**WARNING**: when using hybrid functional, one must copy **OUT.HSEWR($i$)** files from previous SCF calculation to the current NONSCF directory.
{: .block-warning }

#### JOB = DOS
Do density of state calculation.

|||
| --- | --- |
|**Related tags**| DOS\_DETAIL, IN.WG, ... |

DOS calculation requires input wave function and eigen energy, usually from a previous SCF or NONSCF calculations. One must set **IN.WG = T** in **etot.input**. one should also copy or link **OUT.EIGEN** file from previous calculation to the current DOS directory.

DOS uses input wave functions to calculate their projections on atomic orbitals, the nonlocal potential projector in this step is different from NONSCF and SCF. One can get partial and projected DOS by using **plot\_DOS\_interp.x**.

There are two ways to calculate DOS, one is conventional DOS calculation, the other is k-point interpolation scheme for DOS calculation, which can get a smooth DOS with very few k-points.

Frequently used `etot.input` settings for conventional DOS calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = DOS
    IN.PSP1 = Si.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 50
    MP_N123 = 9 9 9 0 0 0 #keep it consistent with previous calculation
    IN.WG = T
```

For k-point interpolation scheme method, it is controled by tag **DOS\_DETAIL**.

Frequently used `etot.input` settings for k-point interpolation scheme DOS calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = DOS
    IN.PSP1 = Si.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 50
    MP_N123 = 9 9 9 0 0 0 #keep it consistent with previous calculation
    DOS_DETAIL = 1 9 9 9
    IN.WG = T
```
>
>**WARNING**: one must copy **OUT.EIGEN** file from previous calculation to the current DOS directory.
{: .block-warning }

#### JOB = MOMENT
Do momentum matrix calculation. Calculates the momentum matrix (oscillator strength) between Kohn-Sham orbitals, and the nonlocal potential is considered.

|||
| --- | --- |
|**Related tags**| IN.WG, ... |


MOMENT calculation requires input wave function, usually from a previous SCF or NONSCF calculations. One must set **IN.WG = T** in **etot.input**.

To calculate the optical absorption spectrum, or dielectric constant, the momentum matrix between Kohn-Sham orbitals $\{ \psi_i \}$ needs to be calculated. Formally, this momentum matrix can be expressed as: $M_x(i,j)=<\psi_i\|P_x\|\psi_j>=-<\psi_i\|{\partial H[k]/\partial k_x}\|\psi_j>=-i <\psi_i\|[H,r_x]\|\psi_j>$. here subscript x actually stands for x, y, z directions. So, there are three matrix (in Cartesian coordinates). $P_x$ is the momentum operator. In the case there is no nonlocal potential, $P_x=i \nabla_x$. If only the $i \nabla_x$ is needed in the calculation, one can use utility **ug\_moment.x** to calculate the matrix based on the output wave function **OUT.WG**.

However, if the nonlocal potential needs to be taken into account, there is an additional term $i(V_{NL}r_x - r_x V_{NL})$, which cannot be calculated easily. The **JOB=MOMENT** is to solve this problem, to include this additional term. The resulting momentum matrix is output in **OUT.MOMENT\_EXT\_KPT**. For example, it can be used for RPA calculation for absorption spectrum or dielectric constant calculations. Including this nonlocal term can increase the oscillator strength $\|M_x\|^2$ by  about $10\%$.

The calculated momentum matrix in **JOB = MOMENT** will be stored in output file **OUT.momentK.($x$).1**.

Frequently used `etot.input` settings for MOMENT calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = MOMENT
    IN.PSP1 = Si.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 50
    MP_N123 = 9 9 9 0 0 0 #keep it consistent with previous calculation
    IN.WG = T
```

#### JOB = RELAX
Do atomic position relaxations and cell relaxation using DFT force and total energy.

|||
| --- | --- |
|**Related tags**| RELAX\_DETAIL, RELAX\_HSE, ... |

Inside the atom.config, in the POSITION section, the last three columns 1,1,1, determine whether this atom will move in the x,y,z directions (see section [atom.config](#atomconfig)): 1,1,1, means move, 0,0,0 means fix. Similarly, if cell relaxation is specified in **RELAX\_DETAIL**, then a **STRESS\_MASK** section in **atom.config** can be used to specified whether one wants to relax all components of the unit cell vector, or only selective components of the cell. Also, if cell relaxation is specified, during the relaxation, the number of plane wave G-vectors are kept unchanged. As a result, after the relaxation, if one wants to redo a calculation with the same Ecut, then due to the change of the cell, the number of G-vector will be different, and the energy, stress etc might be different from the previous relaxation runs. So, one might want to do a relaxation again. Or, one can use **STRESS\_CORR** (see section [STRESS\_CORR](#subsection:CORR)) to make a correction for stress calculation, taking into account the effect of finite Ecut to the calculation of stress.

Each atomic relaxation step will do one SCF calculation. Optionally you can also have **RELAX\_DETAIL** (for general RELAX) and **RELAX\_HSE** (for RELAX in the case of HSE calculation) in the **etot.input**. See below [RELAX\_DETAIL](#section:relaxdetail).

A concise result will be reported in **RELAXSTEPS**. The atomic movements for each relaxation step will be reported in **MOVEMENT**, and final atomic configuration is reported in **final.config**.

Frequently used `etot.input` settings for atomic relaxation calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = RELAX
    RELAX_DETAIL = 1 500 0.01
    IN.PSP1 = Si.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 50
    Ecut2 = 200
    MP_N123 = 9 9 9 0 0 0 #modify ``NK1 NK2 NK3''  according to structure lattice
```

In addition, some specific RELAX parameters can be set. Please refer to section [in.relaxopt](#otherinput:in.relaxopt) for details.

Frequently used `etot.input` settings for cell relaxation calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = RELAX
    RELAX_DETAIL = 1 500 0.01 1 0.01
    IN.PSP1 = Si.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 70
    Ecut2 = 280
    MP_N123 = 9 9 9 0 0 0 #modify ``NK1 NK2 NK3''  according to structure lattice
```

#### JOB = EGGFIT
This is a way to fix the ``egghead'' problem in atomic relaxation.

|||
| --- | --- |
|**Related tags**| EGG\_DETAIL, EGG\_CORR |

The egghead problem is caused by the numerical discretization of the real space using grid **n1 * n2 * n3**. As a result, the atom can have an artificial force towards or away from the grid point. In many cases this problem can cause the system relaxing very slowly when the force is small, or the relaxation energy curve become not smooth. In most cases, this problem can be removed by using **Ecut2=4Ecut**, **Ecut2L=Ecut2** (for norm conserving pseudopotential, NC-PSP). However, sometime even this cannot remove the ``egghead'' problem. In that case, the **Ecut2=4Ecut**, **Ecut2L=4Ecut2** will almost always remove the egghead problem. But the ``Ecut2=4Ecut,Ecut2L=4Ecut2'' could be rather expensive. To keep the calculation in ``Ecut2=4Ecut, Ecut2L=Ecut2'' (for NC-PSP), we provide a **JOB=EGGFIT** procedure to remove the egghead problem. This can be useful for large system relaxation runs. In order to use this procedure, one needs to do the relaxation in two steps:

1. set the **JOB = EGGFIT** in the **etot.input** and give an additional setting: **egg\_detail = np1, np2, np3**; **ECUT2 = 4ECUT**, **ECUT2L = ECUT2**; Here, **np1**, **np2**, **np3** indicates the point to probe inside a grid, usually they are 2,2,2 or 4,4,4. After running PWmat, it will give a new file **CC.egghead** which will be used in the following step.

    Frequently used `etot.input` settings for eggfit calculation:

    ```bash
        4 1
        IN.ATOM = atom.config
        JOB = EGGFIT
        EGG_DETAIL = 2 2 2
        IN.PSP1 = Si.SG15.PBE.UPF
        XCFUNCTIONAL = PBE
        Ecut = 50
        Ecut2 = 200
        MP_N123 = 9 9 9 0 0 0 #modify ``NK1 NK2 NK3''  according to structure lattice
    ```

2. set the **JOB = RELAX** with an additional setting: **EGG\_CORR = T**, **ECUT2 = 4ECUT**, **ECUT2L = ECUT2**. **EGG\_CORR = T** means PWmat will read **CC.egghead** to do the egghead correction during relaxation.

    Frequently used `etot.input` settings for atomic relaxation with egg\_corr:

    ```bash
        4 1
        IN.ATOM = atom.config
        JOB = RELAX
        RELAX_DETAIL = 1 500 0.01
        IN.PSP1 = Si.SG15.PBE.UPF
        XCFUNCTIONAL = PBE
        Ecut = 50
        Ecut2 = 200
        EGG_CORR = T  #read ``CC.egghead'' file from previous ``JOB=EGGFIT''
        MP_N123 = 9 9 9 0 0 0 #modify ``NK1 NK2 NK3''  according to structure lattice
    ```

#### JOB = MD
Do Born-Oppenheimer molecular dynamics (MD) simulations.

|||
| --- | --- |
|**Related tags**| MD\_DETAILS, IN.MDOPT, .. |

Must have variable ``MD\_DETAIL'' in etot.input [MDDETAIL]().
The PWmat can perform: Verlet, Nose-Hoover, Langevin, Berendsen dynsmics. We have a concise output as reported in MDSTEPS. The atomic movements for every step are reported in MOVEMENT. One can also do special MD, e.g., with applied different force on each atom, or different specified temperature on each atom within the Langevin dynamics. Inside the atom.config, in the POSITION section, the last three column 1,1,1, determines whether this atom will move in the x,y,z directions (see section \ref{inputfile:atomconfig}): 1,1,1, means move, 0,0,0 means fix. If dynamics can change the unit cell vector (e.g., in NPT calculation), one can also use STRESS\_MASK section in atom.config to specify which cell vector component to change.

Frequently used `etot.input` settings for MD calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = MD
    MD_DETAIL = 1 1000 1 300 300
    IN.PSP1 = Si.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 50
    MP_N123 = 1 1 1 0 0 0 2
```

In addition, some specific MD parameters can be set. Please refer to section [in.mdopt](#otherinput:in.mdopt) for details.

#### JOB = TDDFT
Do real-time time-dependent DFT calculation (rt-TDDFT).

|||
| --- | --- |
|**Related tags**| MD\_DETAIL, TDDFT\_DETAIL, TDDFT\_TIME, TDDFT\_SPACE, IN.A\_FIELD, TDDFT\_BOLTZMANN, IN.TDDFTOPT, .. |

This is a major functionality in PWmat. It uses a new algorithm as reported in Ref.\cite{pwmat3}. The detail of this **JOB** is described separately in [appendix B](#appendix:tddft). The rt-TDDFT can be used to simulate the dynamic process where both nuclei and electron movements are important, and the electron is no longer in the ground state during the nuclear movement (for example, in a high speed ion collision with a material). It can also be used to study optics (absorption spectrum, or nonlinear optics). It includes both the electron-electron interaction, and electron-phonon interaction. The TDDFT simulation is more expensive than the NAMD calculation. Mostly this is because one needs to use a smaller time step **dt** (e.g., 0.1 fs), and calculate more electron adiabatic states (to expand the time evolving wave functions).

#### JOB = NAMD
Do non-adiabatic molecular dynamics.

|||
| --- | --- |
|**Related tags**| MD\_DETAIL, NAMD\_DETAIL, TDDFT\_TIME, TDDFT\_SPACE, TDDFT\_STIME, IN.A\_FIELD, IN.MDOPT, .. |

This is done under the approximation of Born-Oppenbeimer MD (BO-MD) for nuclear movement. So the time to do NAMD is almost the same as that for MD. The actual NAMD simulation is done as a post-processing after the DFT BO-MD. It will generate a file OUT.NAMD. Some post-processing program (e.g., the "namd\_dm.x" in [Boltzman-NAMD](http://www.pwmat.com/module-download) can be used to study the single carrier dynamics during the BO-MD process. It only simulates the behavior of a single carrier. While it takes into account the effects from other electron and phonon to the dynamics of this single carrier (hence, include the electron-phonon coupling etc), it ignores the effects of this carrier to the dynamics of other electron and phonon (i.e, there is no feedback from carrier to phonon, or carrier to other electron, thus it cannot be used to study polaron effect).
    
Advantageously, it also does not have the erroneous carrier self-interaction. It is suitable to study the carrier dynamics (e.g., charge transfer between molecule, or spin dynamics of one defect) of some large systems. Compare to TDDFT, one advantage is that it can do much bigger system with much longer time. The details are also described in [NAMD](appendix B) setction.

#### JOB = NEB

Do nudged elastic band (NEB) calculation.

|||
| --- | --- |
|**Related tags**| NEB\_DETAIL, .. |

It must have a variable **NEB\_DETAIL** in **etot.input** (see section [NEB\_DETAIL](#subsection:NEBDETAIL)). Besides **IN.ATOM**, which gives the first valley site atomic position, there must be a second valley site position given in the **NEB\_DETAIL** line. One must precalculate (e.g., using **JOB=RELAX**) the atomic configuration of these two valley sites before using **JOB=NEB** to calculate their barrier. See [NEB\_DETAIL](#subsection:NEBDETAIL) for more details and how to set up the calculations. Output files: **RELAXSTEPS**, **NEB.BARRIER**, **MOVEMENT**. **NEB.BARRIER** gives the barrier height information, while **MOVEMENT** gives all the image atom.config files within each NEB step.

>
>**WARNING**: During NEB calculation, if you encounter the following error: ``equivalent atom not found under symm op, stop'', please turn off the symmetry, just set **MP\_N123 = NK1 NK2 NK3 0 0 0 2** (see section [MP\_N123](#subsection:MPN123))
{: .block-warning }

An example etot.input for NEB calculation:

```bash
    4 1
    IN.ATOM = atom1.config
    JOB = NEB
    NEB_DETAIL = 5 100 0.03 5 1 2 -7946.015 -7946.015 1 atom2.config
    ACCURACY = High
    IN.PSP1 = C.SG15.PBE.UPF
    IN.PSP2 = Li.SG15.PBE.UPF
    XCFUNCTIONAL = PBE
    Ecut = 50
    MP_N123 = 1 1 1 0 0 0 2
```

Additional relaxation settings can be found in section [in.relaxopt](#otherinput:in.relaxopt).

#### JOB = DIMER

Do dimer method calculation.

|||
| --- | --- |
|**Related tags**| IN.RELAXOPT, DIMER\_DIR\_N, .. |

Dimer method is used for finding saddle points without knowledge of the final state of the transition is described, and allows users to search for a nearby saddle point from a given initial configuration. The dimer method is designed to deal with problems with unkown reaction mechanisms.

Some specific DIMMER parameters can be set in file **IN.RELAXOPT**, please refer to section [in.relaxopt](#otherinput:in.relaxopt) for details.

The initial direction along the dimer can be set in structure file by tag **DIMER\_DIR\_N** (see section [atom.config](#atomconfig)).

The final configuration is writen in file **final.config**, the configurations of each translation step are in file **MOVEMENT**.

Another file **DIMERSTEPS** can be used to check the convergency, see section [DIMERSTEPS](#outputfile:DIMERSTEPS).

An example etot.input for DIMMER calculation:

```bash
    4 1
    IN.ATOM = atom.config
    JOB = DIMER
    IN.PSP1 = H.SG15.PBE.UPF
    IN.PSP2 = N.SG15.PBE.UPF
    Ecut = 50
    Ecut2 = 200
    MP_N123 = 2 2 2 0 0 0 2
    #in somecases you should turn off the symmetry
    fermidE = 0.2
```

#### JOB = SCFEP

Do electron-phonon coupling calculation.

|||
| --- | --- |
|**Related tags**| SCFEP\_DETAIL, .. |


The procedure is the following, one first carries out a **JOB = SCF** calculation, and **OUT.FORCE = T**, so there will be atomic forces (or perhaps before that, there will be a **JOB = RELAX**, to relax the atoms). Copy **OUT.FORCE** file into **IN.FORCE** (and set **IN.FORCE = T**). For **JOB = SCFEP**, we must have **IN.WG = T**, here **IN.WG** is copied from the **OUT.WG** from the previous **JOB = SCF** calculation. Now, in **JOB = SCFEP** calculation, the state $\psi(ist1,ikpt,ispin)$, and $\psi(ist2,kpt,ispin)$ (to be specied in the line of **SCFEP\_DETAIL**) will be used to calculate: $<\psi(ist1,ikpt,ispin)\|\delta{H}/\delta{R}\|\psi(ist2,ikpt,ispin)>$, this result will be represented as an perturbed atomic forces (the perturbation is proportional to $\alpha$ as specified in **SCFEP\_DETAIL**), and reported in **OUT.FORCE**. The actual coupling constants will be reported in **OUT.EP\_COEFF** (here the **FORCE\_new** has already be subtracted by **IN.FORCE**, and divided by $\alpha$).

The $\alpha$ should be small, something like 0.1, 0.2.

The electron-phonon coupling constant reported in **OUT.EP\_COEFF**, together with phonon calculations can be used to study non-adiabatic decay and charge trapping by defect states.


#### JOB = POTENTIAL

Do potential calculation.

|||
| --- | --- |
|**Related tags**| IN.RHO, .. |

This take the input charge, **IN.RHO = T**, output the potential: **OUT.VR**, **OUT.VR\_HION**, then stop. It will be useful for charge patching and defect calculation. Basically it is a simple Poisson solver. In particular, for isolated systems, one can specify **COULOMB = 1** for Poisson solver without periodic image potential. This can be a very quick calculation.

#### JOB = HPSI

Do H$\psi$ calculation.

|||
| --- | --- |
|**Related tags**| IN.WG, IN.VR, IN.RHO, .. |

This is used to calculate $H\psi_i$ and output the wave function $H\psi_i$ in **OUT.HPSI** (and **OUT.HPSI\_2** for spin=2). It must has a **IN.WG = T**, and have **IN.VR = T** or **IN.RHO = T** to have the proper Hamiltonian. This is provided, so one can carry out some analysis, for example to calculate the electron-phonon coupling.

#### JOB = WKM

Do Wannier Koopmann method (WKM) calculation.

|||
| --- | --- |
|**Related tags**| IN.WANNIER, IN.S\_WKM, .. |


When we do a DFT calculations, one difficulty is to calculate a right bandstructure with right band gaps which can agree well with experimental results. However, it is a common sense that LDA, PBE and even HSE (when $\alpha$ is set to be 0.25) ofter nderestimate band gap results. For band gap calculation, we often define the band gap as the difference between the electron affinity (EA) energy and the ionization energy (IE). Here, EA can be expressed as E(N+1)-E(N) and IE as E(N)-E(N-1). N is the number of electrons in the neutral system, and N+1 and N-1 indicate the system has one more or one less electron respectively. E(N) is the self-consistent energy of the system with N electrons.

[WKM.png]

total energy profile describes an LDA total energy calculation and an ``exact'' energy result for an open system

One possible way to overcome this underestimation is to perform a Koopmans condition on normal DFT calculations.

If we plot E(N) as a function of N in LDA calculations, a black parabola in the above figure  can be plotted. However, this seems to be unphysical. If we want to calculate E(N+s), which 0 < s <1 , the total energy of a system with a fractional number of electrons can be defined as a statistical mixture of the N electron and N+/-1 electron state.This leads to a linear segment total energy function of s, which is plotted in red ``exact'' straight lines in the above figure.This linear segment property is also called Koopmans condition.

However, if we just add s electrons in a unit cell, as a result, accordingto Janak’s theory, the total energy difference is the same asKohn-Sham orbital eigen energy. To overcome this problem, we added an electron into a localized Wannier function instead of the extended Kohn-Sham orbitals. As a result, the WKM total energy can beexpressedas:

$$
        E_{WKM}(\{s_k\}) = E_{LDA}(\{s_k\}) + \sum\limits_{k}E_k(s_k)
$$

Here, w indicates the wannier functions. $s_k$ (0 < $s_k$ < 1) indicates the occupation number of this Wannier function. During the LDA calculation, we add or remove the electron from $\phi_k$ of one spin channel and all the other orbitals in this spin channel should be orthogonal to this Wannier function $\phi_k$. All the other orbitals (except this oneWannier function) are variationally changed to minimize the total energy, which results in the ground state energy E$_{LDA}$({$s_k$}). Thus, a simple analytical expression of $E_{LDA}$({$s_k$}) can be writen as

$$
        E_k(s_k) = \lambda_ks_k(1 - s_k)
$$

The $\lambda_k$ can be determined from $E_{LDA}$({$s_k$}) (to make E$_{WKM}$({$s_k$}) a straight line vs. $s_k$). $\lambda_k$ can be calculate by PWmat JOB = WKM mode.

It requires the input of Wannier wave functions (in real space), provided by files in the name of IN.WANNIER\_00001.u, IN.WANNIER\_00002.u IN.WANNIER\_00001.d IN.WANNIER\_00002.d etc. They are written in the same format as the charge density IN.RHO (thus the Wannier wave functions are real). Each of this file contain only one Wannier function $\phi_k$. The information for these Wannier functions are provided in a file called: IN.S\_WKM. They look like

```bash
S_WKM1
2            : the number of Wannier function in up spin
0.5   1.0    : s1_u ss1_u : the occupation of the first up wannier function
0.0   0.0    : s2_u ss2_u : the occupation of the second up wannier function
S_WKM2
1            : the number of Wannier function in down spin
0.9   1.9    : s1_d, ss1_d : the occupation of the first down wannier function
M_FIX_WKM    : This section is optional
nb_fix1,nb_fix2,iflag_wkm_Hxc,nb_exclude_Hxc
```
In the WKM calculation, the Wannier function $\phi_k$ will be occupied according to s1\_u, or s2\_u. The occupation ss1\_u, ss2\_u is used to make the system a full shell. This is only used when there is a special treatment for the exchange-correlatin functional by exclude some core level charge densities. Otherwise, they are not really used.    In the WKM calculation, the other "normal" wave functions $\{ \psi_i \}$ will be orthogonal to the Wannier functions $\phi_k$ included in the IN.S\_WKM file. Their total charge (from $\{ \psi_i \}$) is determined by **NUM\_ELECTRON**. So, the total charge is **NUM\_ELECTRON** plus the s1\_u etc. It is a good idea to always include **NUM\_ELECTRON** in the WKM calculation.

The optional session **M\_FIX\_WKM**  is used for a special exchange-correlation functional treatmenf for the WKM calculation with semicore states (the states very deep in energy).

We found that, in WKM calculation for the lambda for a Wannier function $\phi_k$, it might be necessary to fix some deep level bands (so they do not change during the SCF calculation. These bands are indicated by  [**nb\_fix1**, **nb\_fix2**]$.
    
If **iflag\_wkm\_Hxc=1**, then the bands: **[1, nb\_exclude\_Hxc]** counted from the bottom will not be included in the exchange-correlation function evaluations, and in this case the **ss1\_u** and **ss2\_u** are used to occupy the Wannier function so to get a closed shell structure.

Note, one can also use JOB=WKM to do some other calculations. For example, to fix some wave function without change during SCF, but to relax all the other wave functions, while keeping all of them orthogonal at the same time. It is not straightforward to do due to the real space form in IN.WANNIER\_00001.u etc, but it can be done.

Note, there is another related calculation, that is the SCF WKM calculation. For that calculation, it is done not by JOB=WKM, instead it is done by using JOB=SCF, and XCFUNCTIONAL = LDAWKM, or XCFUNCTIONAL = LDAWKM2.
    
Those are used to carry out SCF WKM calculation when the WKM parameter $\lambda_k$ has already been calculated, or they can be used to carry out linear response WKM calculations. Please check the XCFUNCTIONAL section for that. 

Detailed steps to perform a WKM calculation in PWmat code, please refer to ``\href{http://www.pwmat.com/module-download}{Band structure calculation by WKM}''.

#### JOB = ATOMIC\_ORB

This will calculate the chosen atom's atomic wavefunction specified in its pseudopotential file.
    
Parameter ATOMIC\_ORBITAL\_IATOM\_OUT in etot.input must be set to set the atom index range(the index is from the IN.ATOM configuration file).

The output filenames of atomic wavefunctions in format atomic\_orb\_iatom\_chi\_ichi\_l\_il\_m\_im, iatom is the index of chosen atom, ichi is the index of PP\_CHI in pseudopotential file, il and im are the corresponding quantum numbers.

The output atomic wavefuncions are in the same format of OUT.RHO, one can use convert\_rho\_new.x to convert to xsf format. All datas of atomic wavefunctions are real numbers.

#### JOB = TRANS

Calculate system state $\psi_l$(r) of transport device based on auxiliary function $W_l$(r). Please refer to [pwmat\_transport](http://www.pwmat.com/module-download) for more details.

An example file etot.input:

```bash
4 1
job = trans
in.atom = system.config
in.vr=T 
SCF_ITER0_1 = 1 100000 3 0.0 0.2 1 # must have this line, will only do one iteration with many CG steps, so NITER0_1=1; NLINE0=100000 
num_band=35 # number of $W_l$(r) 
N123=480 96 32 
Ecut=50 
Ecut2=100 
precision=double
wg_error=1.d-5
flag_cylinder=1  # will not use Ecut value for cutoff energy in direction x\\
in.kpt=T 
IN.PSP1 = Cu.FHI.LDA.UPF
IN.PSP2 = S.FHI.LDA.UPF
IN.PSP3 = C.FHI.LDA.UPF
IN.PSP4 = H.FHI.LDA.UPF
```

### ACCURACY

|Tag|ACCURACY|
| --- | --- |
|**Format**|ACCURACY = NORM / HIGH / VERYHIGH|
|**Default**|NORM|

We have introduced three control flags: Accuracy, precision and convergence for the user to easily control different aspects of the calculation. The values of these flags will change the settings of other more detailed parameters. However, one can also set those parameters directly. Those detailed parameters should have higher priority  (if they are explicitly set) than these three control flags.

Control the calculation accuracy, helping to set up the default values for other parameters in etot.input. This parameter will influence the setting of default ECUT/ECUT2 and P123 (for HSE) (see the following).

```bash
ACCURACY      NORM            HIGH         VERYHIGH
ECUT          PSP/INPUT       PSP/INPUT    PSP/INPUT
ECUT2         2*ECUT          4*ECUT       4*ECUT
ECUT2L        ECUT2(NCPP)     ECUT2        4*ECUT2
              4*ECUT2(USPP)   4*ECUT2      4*ECUT2
ECUTP         ECUT            4*ECUT       4*ECUT
RCUT          PSP/INPUT       PSP/INPUT    1.1*PSP/INPUT
```

**ACCURACY = NORM**, the default ECUT will be used, and ECUT2 = 2 * ECUT, ECUT2L = ECUT2 for NCPP, and ECUT2L = 4 * ECUT2 for ultrasoft PSP. P123 = NP1, NP2, NP3, which equals 2/3 of N1, N2, N3 (or generated from ECUTP, using a FFT box just containing the ECUTP sphere). 

**ACCURACY = HIGH**, if ECUT/ECUT2 are not specified, it will set ECUT = 1.0 * default value in the pseudopotential file and ECUT2 = 4 * ECUT, ECUT2L = ECUT2, and P123 = N123. 

**ACCURACY = VERYHIGH**, if ECUT/ECUT2 are not specified, it will set ECUT = 1.0 * default value in the pseudopotential file and ECUT2 = 4 * ECUT, ECUT2L = 4 * ECUT2, and P123 = N123. 


### PRECISION

|Tag|PRECISION|
| --- | --- |
|**Format**|PRECISION = AUTO / SINGLE / DOUBLE / MIX|
|**Default**|AUTO|

The precision controlling flag of GPU calculation.
    
```bash
PRECISION             AUTO(DEFAULT)   DOUBLE         SINGLE         MIX
SCF(HSE)              NCPP:DOUBLE     NCPP:DOUBLE    NCPP:DOUBLE    NCPP:DOUBLE
                      USPP:SINGLE     USPP:SINGLE    USPP:SINGLE    USPP:SINGLE
RELAX_HSE(NUM_LDA>0)  LDA:SINGLE      LDA:DOUBLE     LDA:SINGLE     LDA:MIX
                      HSE:DOUBLE      HSE:DOUBLE     HSE:DOUBLE     HSE:DOUBLE
SCF,RELAX(LDA/GGA)    SINGLE          DOUBLE         SINGLE         MIX
```
                
**PRECISIION = AUTO**, double or single precision in the calculation will be automatically adjusted.

**PRECISION = SINGLE**, use single precision of GPU calculation, default (except for HSE). For most cases, SINGLE is good enough. Typically, it can converge the total energy to 0.1 meV, and the error for rho to be about 1.E-5, and error for total energy to be about 1.E-4 (eV).

**PRECISION = DOUBLE**, use double precision of GPU calculation. Default for the PBE part of HSE calculation. This however can be slower on the Mstation. Use this only you really want to make sure the numerical precision is not a problem.

**PRECISION = MIX**, use both double and single precisions in the calculation,automatically adjust. Only some critical parts use DOUBLE, the other parts use SINGLE. It is a compromise between SINGLE and DOUBLE precisions. Usually this should be good enough for almost any calculations.

Obviously, from **SINGLE**, **MIX** to **DOUBLE**, more accurate, but more costly. In most calculations, SINGLE is good enough, and MIX can be almost as good as the DOUBLE precision. If there are some issues in terms of convergence and the final result, one can use DOUBLE precision to check. Note, for HSE, the HSE Fock exchange term is calculated with SINGLE, but the SCF iterations are done using DOUBLE.

### CONVERGENCE

|Tag|CONVERGENCE|
| --- | --- |
|**Format**|CONVERGENCE = EASY / DIFFICULT|
|**Default**|EASY|
|**Related tags**|RHO\_RELATIVE\_ERROR, WG\_ERROR, RHO\_ERROR, SCF\_ITER0, SCF\_ITER1, ACCURACY|

Control the convergence parameters of the SCF self-consistent iteration.
    
```bash
CONVERGENCE            EASY                DIFFICULT
WG_ERROR               1.0E-4              0.5E-4
E_ERROR                1.0E-7*TOTNEL*Har   0.01E-7*TOTNEL*Har
RHO_ERROR              0.5E-4              0.5*0.5E-4
RHO_RELATIVE_ERROR     0.0                 0.0
(TOTNEL: total number of electrons)
(Har: 27.21138602eV)
```        

**CONVERGENCE = EASY**, use less self-consistent iteration steps to do the calculation in default setting.%: e.g.: %SCF\_ITER0 = 40 and SCF\_ITER1 = 40. For the normal calculation, we recommend to use this setting. In some cases, it is hard to make the self-consistent iteration converge, you can try the ``DIFFICULT'' value. 

**CONVERGENCE = DIFFICULT**, decrease several parameters for a better SCF convergence.

### NUM\_MPI\_PER\_GPU

|Tag|NUM\_MPI\_PER\_GPU|
| --- | --- |
|**Format**|NUM\_MPI\_PER\_GPU = N|
|**Default**|NUM\_MPI\_PER\_GPU = 1|
This parameter is used to control how many threads are bound to a GPU at the same time.

### NUM\_BLOCKED\_PSI

|Tag|NUM\_BLOCKED\_PSI|
| --- | --- |
|**Format**|NUM\_BLOCKED\_PSI = T/F|
|**Default**|NUM\_BLOCKED\_PSI = F|

In NUM\_BLOCKED\_PSI = T, PWmat will divide the wavefunctions into N parts and then put the parts into GPU memory successively one after another during scf iteration. This is to save the use of GPU memory. If a previous run found the GPU out of memory, this can be tried.

>
>**WARNING** : It is notallowed to use this parameter when ``ENERGY\_DECOMP = T'' in etot.input. If NUM\_BLOCKED\_PSI = T,the decomposed energy can suddenly be very wrong.
{: .block-warning }

This parameter intends to save the GPU memory to calculate a larger or more complicated systems. So when PWmat tells {\color{blue}{``CUDA MEMORY INSUFFICENT''}}, one can try this parameter by setting NUM\_BLOCKED\_PSI=T, and etc. Note that using this parameter will reduce the speed of PWmat (e.g., by a factor of 1.5).

### WF\_STORE2DISK

|Tag|WF\_STORE2DISK|
| --- | --- |
|**Format**|WF\_STORE2DISK = 1 / 0|
|**Default**|WF\_STORE2DISK = 0|

If WF\_STORE2DISK = 1, the wavefunctions will be written into disk, otherwise written into cpu memory.
This parameter is used to save cpu memory to calculate a larger or more complicated systems, in particular for the case multiple k-points are calculated (then one can use WF\_STORE2DISK=1). Note that: it will reduce the performance of PWmat in some degree.

### USE\_GAUSSIAN

|Tag|USE\_GAUSSIAN|
| --- | --- |
|**Format**|USE\_GAUSSIAN = T/F EPS\_GAUSSIAN NEIGH\_RADIUS IS\_PERIODIC\_A IS\_PERIODIC\_B IS\_PERIODIC\_C ELPA2\_OR\_1 USE\_PERTURB DN\_PERTURB START\_ITER\_PERTURB|
|**Default**|USE\_GAUSSIAN = F 1.E-8 15.0 F F F 2 F 5 5|

If USE\_GAUSSIAN = T, PWmat will use gaussian basis instead of plane wave basis. EPS\_GAUSSIAN determines the accuracy of the gaussian basis on grid of real space, usually shoule be around 1.E-8 ~ 1.E-10. NEIGH\_RADIUS defines whether two atoms are neighbors, if the distance of two atoms is within the NEIGH\_RADIUS, then these two atoms have interactions, this can be used to reduce the amount of calculations because of the localization of gaussian basis, the unit of NEIGH\_RADIUS is bohr. IS\_PERIODIC\_A, IS\_PERIODIC\_B, IS\_PERIODIC\_C can be T/F, these parameters define whether the structrure is periodic or not along each direction of lattice vector. ELPA2\_OR\_1 can be 2 or 1, if 2, use elpa 2-stage solver; if 1, use elpa 1-stage solver. USE\_PERTURB can be T/F, if T, in the process of SCF, from the SCF iteration of START\_ITER\_PERTURB, when mod(SCF\_iteration, DN\_PERTURB) = 0, or say SCF\_iteration can be divided by DN\_PERTURB, use direct diagonalization to solve the eigen energys and wavefunctions, otherwise use perturbation method.

If USE\_GAUSSIAN = T, there need another two files IN.GAUSSIAN and IN.POTENTIAL. IN.GAUSSIAN specifies the basis of each type of atoms, the format of IN.GASSIAN:
```bash 
        2
        H SZV-GTH-q1 GTH_BASIS_SETS
        C SZV-GTH-q4 GTH_BASIS_SETS
```
The first line is the number of atom types, the lines followed are each type's basis. Start from the second line, the first column is the name of element, the second column is the name of basis, the third column is the name of basis file. The basis file can be downloaded from PWmat's website or other opensource gaussian basis. The order of element type should be the same with IN.PSP* in etot.input.

IN.POTENTIAL specifies the pseudopotentials of each type of atoms, the format of IN.POTENTIAL:
```bash
        2
        H GTH-PBE-q1 GTH_POTENTIALS
        C GTH-PBE-q4 GTH_POTENTIALS
```
The first line is the number of atom types, the lines followed are each type's psedupotential. Start from the second line, the first column is the name of element, the second column is the name of pseudopoential, the third column is the name of psedupotential file. The basis file can be downloaded from PWmat's website or other opensource GTH psedupotentials. When USE\_GAUSSIAN = T, PWmat use both gaussian basis and analytical GTH psedupotentials to do the calculation, but in the current implementation, PWmat still need to set the IN.PSP* which specify the psedupotentials used in plane wave basis implementation, usually the SG15 norm conserving pseudopoentials.

In current version PWmat with gaussian basis can just do a limited JOB=SCF/MD/RELAX, without stresses, HSE, LDA+U, SOC and many other parameters. But you can check the total energy, band structure, charge density and forces, and use multiple k-points.
    
For some structrures with periodic boundaries with small lattices, if you encounter some problems with the eigen energy solving step, you can try to change ELPA2\_OR\_1, this may help, or you can try to use different gaussian basis.

An example of gaussian basis,

file atom.config:

```bash
        5
        Lattice vector
            15.0000000000     0.0000000000      0.0000000000
            0.0000000000     15.0000000000      0.0000000000
            0.0000000000      0.0000000000     15.0000000000
        Position, move_x, move_y, move_z
        6    0.500000000000    0.500000000000    0.500000000000 0 0 0
        1    0.474669993000    0.430790007000    0.518549979000 0 0 0
        1    0.474669993000    0.518549979000    0.430790007000 0 0 0
        1    0.474669993000    0.550670028000    0.550670028000 0 0 0
        1    0.575999975000    0.500000000000    0.500000000000 0 0 0
```

file etot.input:
```bash
        4 1
        job=scf
        in.atom = atom.config
        in.psp1=H.SG15.PBE.UPF
        in.psp2=C.SG15.PBE.UPF
        use_gaussian=T 1.d-10 15.0 F F F 2 F 5 5
```

file IN.GAUSSIAN:
```bash
        2
        H SZV-GTH-q1 GTH_BASIS_SETS
        C SZV-GTH-q4 GTH_BASIS_SETS
```

file IN.POTENTIAL:
```bash
        2
        H GTH-PBE-q1 GTH_POTENTIALS
        C GTH-PBE-q4 GTH_POTENTIALS
```

file GTH\_BASIS\_SETS:
```bash
        ......
        H SZV-GTH-q1 SZV-GTH
        1
        1  0  0  4  1 
        8.3744350009  -0.0283380461
        1.8058681460  -0.1333810052
        0.4852528328  -0.3995676063
        0.1658236932  -0.5531027541
        ......
```

file GTH\_POTENTIALS:
```bash
        ......
        C GTH-PBE-q4 GTH-PBE
        2    2
        0.33847124    2    -8.80367398     1.33921085
        2
        0.30257575    1     9.62248665
        0.29150694    0
        ......
```
Here are some output files when USE\_GAUSSIAN=T: 

|||
| --- | --- |
|**OUT.GAUSSIAN\_H**|The matrix elements of the Hamiltonian in the Gaussian basis.  |
|**OUT.GAUSSIAN\_S**|The matrix elements of the overlap matrix in the Gaussian basis.  |
|**OUT.GAUSSIAN\_H\_T**|The matrix elements of the Hamiltonian in the Gaussian basis, labeled by translation vector T.  |
|**OUT.GAUSSIAN\_S\_T**|The matrix elements of the overlap matrix in the Gaussian basis, labeled by translation vector T.  
|
|**OUT.GAUSSIAN\_BASIS\_INDEX**|The index of the basis function, the position of the basis function, the type of the basis function.|

    
Formula of OUT.GAUSSIAN\_H\_T:

$$H_{m,m'}(T) = \int dr \chi^*_m(r-\tau_m) \hat{H} \chi_{m'}(r-(\tau_{m'}+T))$$, 

where m and m' are the index of the basis, T is the translation vector, $\chi_m$ is the m-th basis function, $\tau_m$ is the position of the m-th basis function(i.e. the position of the atom which the m-th basis function belongs to), $\hat{H}$ is the Hamiltonian operator.You can use following fortran codes to read the file OUT.GAUSSIAN\_H\_T:

```fortran
        subroutine read_gaussian_H_T()
        implicit none
        complex(kind=8), allocatable, dimension(:, :, :, :, :) :: H
        integer :: Nx_pbc_t, Ny_pbc_t, Nz_pbc_t, num_mcgtos_t
        integer :: i, j, k, l, m, n
        real*8 :: AL(3, 3)
        complex(kind=8), allocatable, dimension(:, :) :: T_all
        integer :: T(3)
        !
        open (10, file="OUT.GAUSSIAN_H_T", form="unformatted", status="old")
        read (10) Nx_pbc_t, Ny_pbc_t, Nz_pbc_t, num_mcgtos_t, AL
        allocate (H(num_mcgtos_t, num_mcgtos_t, &
        -Nx_pbc_t:Nx_pbc_t, -Ny_pbc_t:Ny_pbc_t, -Nz_pbc_t:Nz_pbc_t))
        allocate (T_all(num_mcgtos_t, num_mcgtos_t))
        do i = -Nx_pbc_t, Nx_pbc_t
            do j = -Ny_pbc_t, Ny_pbc_t
                do k = -Nz_pbc_t, Nz_pbc_t
                    read (10) T(1), T(2), T(3)
                    do l = 1, num_mcgtos_t
                        read (10) T_all(:, l)
                    end do
                    do l = 1, num_mcgtos_t
                        do m = 1, num_mcgtos_t
                            H(l, m, i, j, k) = T_all(l, m)
                        end do
                    end do
                end do
            end do
        end do
        close (10)
        deallocate (T_all)
        deallocate (H)
        end subroutine read_gaussian_H_T
        !
```

Where AL(1:3,1:3) is the lattice vectors, T(1:3) is the translation vector, num\_mcgtos\_t is the number of the basis functions. (Note, T is an integer array, the acctual translation vector is AL*T.)

Formula of OUT.GAUSSIAN\_S\_T:

$$S_{m,m'}(T) = \int dr \chi^*_m(r-\tau_m) \chi_{m'}(r-(\tau_{m'}+T))$$

where m and m' are the index of the basis, T is the translation vector, $\chi_m$ is the m-th basis function, $\tau_m$ is the position of the m-th basis function(i.e. the position of the atom which the m-th basis function belongs to). You can use following fortran codes to read the file OUT.GAUSSIAN\_S\_T:

```fortran
        subroutine read_gaussian_S_T()
        implicit none
        complex(kind=8), allocatable, dimension(:, :, :, :, :) :: S
        integer :: Nx_pbc_t, Ny_pbc_t, Nz_pbc_t, num_mcgtos_t
        integer :: i, j, k, l, m, n
        real*8 :: AL(3, 3)
        complex(kind=8), allocatable, dimension(:, :) :: T_all
        integer :: T(3)
        !
        open (10, file="OUT.GAUSSIAN_S_T", form="unformatted", status="old")
        read (10) Nx_pbc_t, Ny_pbc_t, Nz_pbc_t, num_mcgtos_t, AL
        allocate (S(num_mcgtos_t, num_mcgtos_t, &
        -Nx_pbc_t:Nx_pbc_t, -Ny_pbc_t:Ny_pbc_t, -Nz_pbc_t:Nz_pbc_t))
        allocate (T_all(num_mcgtos_t, num_mcgtos_t))
        do i = -Nx_pbc_t, Nx_pbc_t
            do j = -Ny_pbc_t, Ny_pbc_t
                do k = -Nz_pbc_t, Nz_pbc_t
                    read (10) T(1), T(2), T(3)
                    do l = 1, num_mcgtos_t
                        read (10) T_all(:, l)
                    end do
                    do l = 1, num_mcgtos_t
                        do m = 1, num_mcgtos_t
                            S(l, m, i, j, k) = T_all(l, m)
                        end do
                    end do
                end do
            end do
        end do
        close (10)
        deallocate (T_all)
        deallocate (S)
        end subroutine read_gaussian_S_T
        !
```

where AL(1:3,1:3) is the lattice vectors, T(1:3) is the translation vector, num\_mcgtos\_t is the number of the basis functions. (Note, T is an integer array, the acctual translation vector is AL*T.)

Formula of OUT.GAUSSIAN\_H:

$$H^\sigma_{m,m'}(k) = \sum_T exp(ik\cdot T) H^\sigma_{m,m'}(T)$$

 where m and m' are the index of the basis, k is the K-point, $\sigma$ is the spin index (Note in current version $\sigma$ can just be 1). You can use following fortran codes to read the file OUT.GAUSSIAN\_H:

```fortran
        subroutine read_gaussian_H()
        implicit none
        complex(kind=8), allocatable, dimension(:, :, :, :) :: H
        integer :: nkpt_t, islda_t, num_mcgtos_t
        integer :: i, j, k, l, m, n, iislda, ikpt,iislda_t, ikpt_t
        real*8 :: AL(3, 3)
        real*8 :: akx_2_t, aky_2_t, akz_2_t
        complex(kind=8), allocatable, dimension(:, :) :: T_all
        !
        open (10, file="OUT.GAUSSIAN_H", form="unformatted", status="old")
        read (10) nkpt_t, islda_t, num_mcgtos_t, AL

        allocate (H(num_mcgtos_t, num_mcgtos_t, nkpt_t, islda_t))

        allocate (T_all(num_mcgtos_t, num_mcgtos_t))
        do iislda = 1, islda_t
            do ikpt = 1, nkpt_t
                read (10) iislda_t, ikpt_t, akx_2_t, aky_2_t, akz_2_t
                do i = 1, num_mcgtos_t
                    read (10) T_all(:, i)
                end do
                do i = 1, num_mcgtos_t
                    do j = 1, num_mcgtos_t
                        H(i, j, ikpt, iislda) = T_all(i, j)
                    end do
                end do
            end do
        end do
        close (10)
        deallocate (T_all)
        deallocate (H)
        end subroutine read_gaussian_H
        !
```
where AL(1:3,1:3) is the lattice vectors, akx\_2\_t, aky\_2\_t, akz\_2\_t are the K-point coordinates, num\_mcgtos\_t is the number of the basis functions. (Note, T is an integer array, the acctual translation vector is AL*T.)

Formula of OUT.GAUSSIAN\_S:

$$S_{m,m'}(k) = \sum_T exp(ik\cdot T) S_{m,m'}(T)$$ 

where m and m' are the index of the basis, k is the K-point.
You can use following fortran codes to read the file OUT.GAUSSIAN\_S:
```fortran
        subroutine read_gaussian_S()
        implicit none
        complex(kind=8), allocatable, dimension(:, :, :, :) :: S
        integer :: nkpt_t, islda_t, num_mcgtos_t
        integer :: i, j, k, l, m, n, iislda, ikpt,iislda_t, ikpt_t
        real*8 :: AL(3, 3)
        real*8 :: akx_2_t, aky_2_t, akz_2_t
        complex(kind=8), allocatable, dimension(:, :) :: T_all
        !
        open (10, file="OUT.GAUSSIAN_S", form="unformatted", status="old")
        read (10) nkpt_t, islda_t, num_mcgtos_t, AL

        allocate (S(num_mcgtos_t, num_mcgtos_t, nkpt_t, islda_t))

        allocate (T_all(num_mcgtos_t, num_mcgtos_t))
        do iislda = 1, islda_t
            do ikpt = 1, nkpt_t
                read (10) iislda_t, ikpt_t, akx_2_t, aky_2_t, akz_2_t
                do i = 1, num_mcgtos_t
                    read (10) T_all(:, i)
                end do
                do i = 1, num_mcgtos_t
                    do j = 1, num_mcgtos_t
                        S(i, j, ikpt, iislda) = T_all(i, j)
                    end do
                end do
            end do
        end do
        close (10)
        deallocate (T_all)
        deallocate (S)
        end subroutine read_gaussian_S
        !
```

Where AL(1:3,1:3) is the lattice vectors, akx\_2\_t, aky\_2\_t, akz\_2\_t are the K-point coordinates, num\_mcgtos\_t is the number of the basis functions. (Note, T is an integer array, the acctual translation vector is AL*T.)

The unit of AL is Bohr, the unit of akx\_2\_t, aky\_2\_t, akz\_2\_t is 2$\pi$/bohr.

    OUT.GAUSSIAN\_BASIS\_INDEX: The index of the basis functions of all atoms, with the following format,

You can use following fortran codes to read the file OUT.GAUSSIAN\_BASIS\_INDEX:

```fortran
        subroutine read_gaussian_basis_index()
        implicit none
        integer :: natom_t, num_mcgtos_t
        integer, allocatable, dimension(:, :) :: basis_range_iatom
        integer :: i, j, k, l, m, n
        !
        open (10, file="OUT.GAUSSIAN_BASIS_INDEX", form="unformatted", & 
              status="old")
        read (10) natom_t, num_mcgtos_t
        allocate (basis_range_iatom(2, natom_t))
        read (10) basis_range_iatom(:, 1:natom_t)
        close (10)
        deallocate(basis_range_iatom)
        end subroutine read_gaussian_basis_index
        !
```

where basis\_range\_iatom(1, i) is the index of the first basis function of the i-th atom, basis\_range\_iatom(2, i) is the index of the last basis function of the i-th atom. Note, the index of atoms may not be the original index of atoms in atom.config(specified by IN.ATOM), PWmat will reorder the atoms, so check output file ORIGIN.INDEX to check index mapping.

## System tags

### ECUT

|Tag|ECUT|
| --- | --- |
|**Format**|ECUT = E|
|**Default**|ECUT = ``WFC\_CUTOFF'' in pseudopotential file|

The plane wave cutoff energy for wavefunction (in $Ryd$, note: $1 Ryd = 13.6057 eV$). The default value of ECUT is taken from the pseudopotential files atom.upf from its WFC\_CUTOFF value. Note, in an plane wave calculation, the plane wave functions with their G-vector ($exp(-iG*x)$) within the energy sphere of ECUT is used as the basis function. Thus, ECUT control the size of the plane wave basis set, is one of the most important calculating parameter.

### ECUT2

|Tag|ECUT2|
| --- | --- |
|**Format**|ECUT2 = E|
|**Default**|ECUT2 = 2*ECUT when ACCURACY = NORM, ECUT2 = 4*ECUT when ACCURACY = HIGH or ACCURACY = VERYHIGH|

The cutoff energy for the soft charge density and the potential (in $Ryd$). In a plane wave calculation, not only the orbital are expanded by the plane waves, the charge density is also expanded by the plane waves. However, the plane wave basis set (within energy ECUT2) used to expand the charge density is larger than the plane wave basis set (within energy ECUT) used to expand the orbital.

Ideally (for high accurate calculations), ECUT2 should equal 4*ECUT. But in reality, smaller ECUT2 can sometime be used, e.g., 3*ECUT, or 2*ECUT. By default, ECUT2 = 2*ECUT for normal accuracy calculation (ACCURACY=NORM), and ECUT2=4*ECUT for high accuracy calculation (ACCURACY=HIGH or VERYHIGH). For JOB=RELAX, to avoid the egghead jittering effect, we recommend to use ECUT2=4*ECUT.

The N1, N2, N3 are determined by ECUT2. Note, the RHO\_CUTOFF value in the pseudopotential files atom.upf **is not used**. Also note that, if N1, N2, N3 are not set (by $N123$), they will be generated by ECUT2, together with NODE1. So, if different number of node NODE1 are used, the N1, N2, N3 values could be different even for the same ECUT2. This is because N1\*N2\*N3 must be evenly divided by NODE1.

### ECUT2L

|Tag|ECUT2L|
| --- | --- |
|**Format**|ECUT2L = E|
|**Default**|ECUT2L = ECUT2 when ACCURACY = NORM or ACCURACY = HIGH, ECUT2L = 4*ECUT2 when ACCURACY = VERYHIGH or use PAW pseudopotential|

The cutoff energy for the hard charge density (in $Ryd$).

Sometime it is necessary to further increase the accuracy of the description for the charge density rho(r) before it is used to calculate the potential via the exchange-correlation functional. Thus, we have a so-called hard charge density, which is described by a plane wave basis set within ECUT2L.

Usually, ECUT2L = ECUT2 for norm conserving pseudopotentials, and ECUT2L = 4 * ECUT2 for ultrasoft pseudopotentials. However, sometime to completely remove the egghead problem, we can also use ECUT2L = 4 * ECUT2 even for the norm conserving psp. Nevertheless, that egghead problem can usually be solved by using EGG\_FIT, so we can still use ECUT2L = ECUT2 for norm conserving pseudopotential (but usually require ECUT2 = 4 * ECUT).

### ECUTP

|Tag|ECUTP|
| --- | --- |
|**Format**|ECUTP = E|
|**Default**|ECUTP = ECUT when ACCURACY = NORM, ECUTP = 4*ECUT when ACCURACY = HIGH or ACCURACY = VERYHIGH|

The cutoff energy to generate P123 for Fock exchange integral evaluation for HSE calculations (in $Ryd$). If P123 is exlicitly input, the P123 will have higher priority.

Note, in order to have accurate results, one needs to have ECUTP = 4*ECUT. This is necessary for accurate force calculations, e.g., during atomic relaxation or phonon mode calculations. Otherwise, the total force might not be zero.

However, if only electronic structure is needed, or for molecular dynamics, or evern TDDFT simulations, one might be able to set a smaller ECUTP, for example, ECUTP=ECUT.

### N123

|Tag|N123|
| --- | --- |
|**Format**|N123 = N1, N2, N3|
|**Default**|N123 = determined by ECUT2|

N1, N2, N3 are the real space grid to describe the wave function or soft charge density in real space. It is also the FFT grid. The default values are determined by ECUT2 (i.e., make sure the ECUT2 sphere can be held inside the $N1, N2, N3$ reciprocal box). Roughly speaking, $N_i=\sqrt{2*ECUT2}/(\pi*\|ALI(:,i)\|)$, here $ALI(:,i)$ is the reciprocal lattice of the input cell lattice $AL(:,i)$ in atom.config file, in $Bohr$ unit, and $\|ALI(:,i)\|$ is the length of the vector. Also, in this formula, $ECUT2$ is in the unit of Hartree.

Note, in our current implementation, $N1 * N2 * N3$ need to be divided evenly by node1. So, sometime it might be necessarily to readjust $N1, N2, N3$ manually (or to change node1). If $N123$ is not explicitly set, their values will be determined automatically by ECUT2 and NODE1. Note, for the same ECUT2, for different NODE1, it can lead to different $N123$.

### N123\_METH

|Tag|N123\_METH|
| --- | --- |
|**Format**|N123\_METH = 0/1|
|**Default**|N123\_METH = 0|

If N123\_METH = 0, check the section N123. If N123\_METH = 1, the way to generate grid of real space will be different with the detault settings, N1\*N2\*N3 will not be required to be divisible by NODE1. Instead, N1\*N2 must be divisible by NODE1 and N1\*(N3/2+1) be divisible by NODE1.


### N123L

|Tag|N123L|
| --- | --- |
|**Format**|N123L = N1L, N2L, N3L|
|**Default**|N123L = determined by ECUT2L|

N1L, N2L, N3L are the real space grid for hard charge density. The default values are determined by ECUT2L. For norm conserving pseudopotential, the soft charge equals hard charge, ECUT2L=ECUT2, so N1L, N2L, N3L equal N1, N2, N3. For ultrasoft, ECUT2L = 4 * ECUT2, N1L, N2L, N3L = 2 * N1, 2 * N2, 2 * N3.

### NS123

|Tag|NS123|
| --- | --- |
|**Format**|NS123 = N1S, N2S, N3S|

N1S, N2S, N3S are the real space grid to calculate the real space nonlocal pseudopotential projector function. So, these are only used for NONLOCAL = 2. For small systems, N1S,N2S,N3S can be larger than N1,N2,N3. For large systems, smaller values can be used to save time for projector generation. Usually these parameters are set automatically.

### MP\_N123

|Tag|MP\_N123|
| --- | --- |
|**Format**|MP\_N123 = NK1, NK2, NK3, SK1, SK2, SK3, FLAG\_SYMM|
|**Default**|MP\_N123 = 1 1 1 0 0 0 0|

This variable is the Monkhorst-Pack grids to generate the reduced k-points. When this line is provided, the PWmat will generate the OUT.SYMM and OUT.KPT using the above Monkhorst-Pack parameters, and the PWmat will continue to run the JOB using these k-points and symmetries.

>
>**WARNING**: if one wants to generate only gamma point, but does not want to use symmetry, one needs to set: ``MP\_N123 = 1 1 1 0 0 0 2''.
{: .block-warning}

Note: if file `IN.KPT` exists and IN.KPT=T, but at the same time, MP\_N123 is also specified, PWmat will ignore **MP\_N123** and use kpoints readin from file `IN.KPT` (i.e, IN.KPT=T has higher priority than MP\_N123). If you want to use the input kpoints and symmetry, you should set IN.KPT=T, IN.SYMM=T, then PWmat will read the files `IN.KPT`, `IN.SYMM` for the calculation.

The SK1, SK2 and SK3 must be either 0 (no offset) or 1 (grid displaced by half a grid point in the corresponding direction). This is the standard options to generate the Monkhorst-Pack k-point grid.

The FLAG\_SYMM controls the symmetry operation(the operations are stored in OUT.SYMM) of k-points. One can refer to the OUT.SYMM for the specific symmetry operations.

FLAG\_SYMM=0, generate kpoints with spatial symmetry and time reversal symmetry. This flag will have the full symmetry operations.

FLAG\_SYMM=1, generate kpoints with spatial symmetry but no time reversal symmetry. This may generate lower symmetry than flag=0. This for example, can be used for magnetic system calculation and for systems with an external magnetic field.

FLAG\_SYMM=2, generate kpoints without any symmetry, i.e. the symmetry is identity operation. This for example can be used for rt-TDDFT simulation with external potential.

FLAG\_SYMM=3, generate kpoints with time reversal symmetry but no spatial symmetry. This may generate lower symmetry than flag=0. This for example, can be used for rt-TDDFT without magnetic moment.

Special attention needs to be paid for the symmetry operation when magnetic system is calculated, e.g., when calculating antiferromagnetic system, since the symmetry operation might ignore the magnetic moment difference between different atoms.

In above, the symmetry includes both point group symmetry operations and space group symmetry operations. Also, for best point group symmetry, one should always place the high symmetry point at the origin (0,0,0) position, since that is the symmetry operation point.

There are several issues one needs to know. If IN.VEXT is used, it is the user's responsibility to figure out what symmetry one should use, since when the above symmetry operation are generated, it does not consider IN.VEXT. For system with MAGNETIC moment section in atom.config and SPIN=2, the MAGNETIC moment of each atom is also used to figure out the symmetry (e.g., one atom with magnetic moment 2 will be different from another atom with same atomic number, but magnetic moment equals -2).

Notes about IN.KPT, IN.SYMM and MP\_N123:

1. Default: IN.KPT=F, IN.SYMM=F, MP\_N123= 1 1 1 0 0 0; IN.KPT=T has higher priority then MP\_N123;

2. IN.KPT=T: read kpoints from file IN.KPT; write OUT.KPT;

3. IN.KPT=F: use MP\_N123 to generate kpoints; write OUT.KPT;

4. IN.SYMM=T and IN.KPT=T: read symmetry operations from file IN.SYMM; write OUT.SYMM;

5. IN.SYMM=T and IN.KPT=F: read symmetry operations from file IN.SYMM, then MP\_N123 will use these symmetry operations to generate kpoints; write OUT.SYMM;

6. IN.SYMM=F and IN.KPT=T: no symmetry used; write 'identity' operation to file OUT.SYMM;

7. IN.SYMM=F and IN.KPT=F: use MP\_N123 to generate kpoints and symmetry operations; write OUT.KPT and OUT.SYMM;

8. for JOB={MD,NEB,TDDFT,NAMD} or SPIN=222: default FLAG\_SYMM = 2.

### SYMM\_PREC

|Tag|SYMM\_PREC|
| --- | --- |
|**Format**|SYMM\_PREC = distance\_tolerance|
|**Default**|SYMM\_PREC = 1.0E-5|

The distance tolerance for symmetry operation. If the distance between two atoms is smaller than SYMM\_PREC, they are considered as the same atom. The default value is 1.0E-5.

### P123

|Tag|P123|
| --- | --- |
|**Format**|P123 = NP1, NP2, NP3|
|**Default**|generated by ECUTP|

When using HSE method, a small box FFT (with grid NP1,NP2,NP3) can be used to calculate the explicit FOCK exchange integral if only electronic structure is needed. This can significantly speedup the calculation without much loss of accuracy (for electronic structure, molecular dynamics, or TDDFT). Sometime NP1, NP2, NP3 can be as small as half of N1, N2, N3. NP1,2,3 are generated using ECUTP. However, in order to have higher accuracy force and energy, one should used ECUTP=4*ECUT, or usually at least ECUTP=ECUT2 (e.g., for atomic relaxation or phonon mode calculations).

### SPIN

|Tag|SPIN|
| --- | --- |
|**Format**|SPIN = 1 / 2 / 22 / 222|
|**Default**|SPIN=1|

SPIN = 1, non-spin-polarized calculation (default). Each orbital will be occupied by 2 electron.

SPIN = 2, spin-polarized calculation, LSDA (magnetization along z axis). For systems except ferromagnetic system, please specify the initial magnetic moment in `atom.config` with the tag section: **MAGNETIC**. Note, for SPIN=2, IN/OUT charge density will have: IN.RHO, IN.RHO\_2, and OUT.RHO, OUT.RHO\_2 (spin up and down components). Similarly, IN/OUT potential will also have spin up and down components: IN.VR, IN.VR\_2, OUT.VR, OUT.VR\_2.

SPIN = 22, spin-orbit coupling (SOC) calculation, but without magnetic moment. This is suitable for semiconductors like CdSe. In this case, each orbital will have spin-up and spin-down components (spinor). But there are also spin-up orbital and spin-down orbital (each of them has spin-up and spin-down components), and both are occupied. As a result there is no magnetic moment. Since there is no magnetic moment, the charge density and potential will only have one component, or say the spin up and down components are the same. So, there will be IN.RHO,OUT.RHO,IN.VR,OUT.VR, but there will be no IN/OUT.RHO\_2, IN/OUT.VR\_2 counterparts.

SPIN = 222, spin-orbit coupling calculation, with noncollinear magnetization in generic directions. For SPIN=222, please specify the initial magnetic moment in `atom.config` with the tags **MAGNETIC\_XYZ**. In this case, the IN/OUT.RHO will also have IN/OUT.RHO\_SOM (a complex 2x2 spin matrix density). IN/OUT.VR will also have IN/OUT.VR\_SOM (a complex 2x2 spin matrix potential) and IN/OUT.VR\_DELTA (a real up-down diagonal potential, note, there is no IN/OUT.RHO\_DELTA).

For SPIN=22 and SPIN=222, the SOC pseudopotentials need to be used. Check the pseudopotential sets, choose the proper SOC pseudopotentials. Note, in our calculation, all the SOC comes from the core levels, so only the heavy atoms have the SOC pseudopotentials. We do not calculate the valence band SOC. For light elements like C, N, O, there is no SOC pseudopotential. However, the SOC pseudopotential and non-SOC pseudopotential can be used in mix within a single calculation. 

>
>WARNING: for SPIN=22 or SPIN=222, you must set parameter **ECUT** in etot.input file, and do not use the default value in pseudopotential file.
{: .block-warning}

SOC can also used together with HSE calculations. This is important for topological system calculations.

### NUM\_ELECTRON

|Tag|NUM\_ELECTRON|
| --- | --- |
|**Format**|NUM\_ELECTRON = value|
|**Default**|NUM\_ELECTRON = value for neutral system|

The total number of occupied valence electron in the system. One can use this to make the system charged, or not charged. Note, for charged system calculations, a uniformed back ground charge is used to solve the Possion equation for COULOMB=0. Default value is the value for neutral system.

### NUM\_ELECTRON\_SPIN

|Tag|NUM\_ELECTRON\_SPIN|
| --- | --- |
|**Format**|NUM\_ELECTRON\_SPIN = NUM\_UP  NUM\_DN|

The number of spin-up and spin-down electrons. This is used for spin-polarized calculation (SPIN=2). If NUM\_ELECTRON\_SPIN is explicitly set in etot.input, it will seperately fix the spin-up and spin-down number of electrons to NUM\_UP and NUM\_DN.

### NUM\_BAND

|Tag|NUM\_BAND|
| --- | --- |
|**Format**|NUM\_BAND = value|
|**Default**|NUM\_BAND = min[1.05*NUM\_ELECTRON/2+10] (SPIN = 1) NUM\_BAND = min[1.2*min[1.05*NUM\_ELECTRON/2+10]] (SPIN = 2) NUM\_BAND = min[1.05*NUM\_ELECTRON+10] (SPIN = 22/222)|

 The number of orbitals to be calculated. When SPIN=2, there are NUM\_BAND spin-up orbitals and NUM\_BAND spin-down orbitals.

### RCUT

|Tag|RCUT|
| --- | --- |
|**Format**|RCUT = value|
|**Default**|RCUT = max[IN.PSP\_RCUT1, IN.PSP\_RCUT2, ..., IN.PSP\_RCUTi]|

The RCUT (in $Bohr$ unit, note: $1 Bohr = 0.529177 \times 10^{-10} m$) is for the cut off radius for nonlocal pseudopotential implementations. It defines the core radius of the nonlocal part. If RCUT is not specified, PWmat will use the maximum of IN.PSP\_RCUTi as the value of Rcut.


### IN.PSP\_RCUT

|Tag|IN.PSP\_RCUT|
| --- | --- |
|**Format**|IN.PSP\_RCUTi = value|

The Rcut of each element type. The IN.PSP\_RCUTi are in the unit of Bohr, not Amgstron. This provides a way to selectively choose the Rcut for different atoms for the nonlocal pseudopotentials. **i** is the index of element  type, in the same order of IN.PSP.

If the IN.PSP\_RCUTi are not specified, the value of IN.PSP\_RCUTi is taken from the pseudopotential file from its **Rcut** value. If these IN.PSP\_RCUTi are not provided in the pseudopotential file, the default value will be set to 3.5. Note, to get a smooth force and energy (e.g., for good RELAX convergence), sufficiently large rcuti should be used (e.g., 4.0 Bohr). However, if Ecut2 is very large, then a slightly smaller rcuti can be used.

### SOM\_SPHERE\_RCUT

|Tag|SOM\_SPHERE\_RCUT|
| --- | --- |
|**Format**|SOM\_SPHERE\_RCUT = value|
|**Default**|SOM\_SPHERE\_RCUT = RCUT|

SOM\_SPHERE\_RCUT is used to determine the spin component for each atom. Roughly, it should be half the bond length (in $Angstrom$). 

>WARNING: The default value is large, and you usually need to set a appropriate value according to the bond length.
{: .block-warning}

## Electronic tags

### E\_ERROR

|Tag|E\_ERROR|
| --- | --- |
|**Format**|E\_ERROR = value|
|**Default**|E\_ERROR = 2.7211E-6*(total number of electrons)|

The error tolerance (convergence criterion) for the total energy (Hartree) in the SCF iterations. The default value is: 2.7211E-6*(total number of electrons) ($eV$). This is related to the SCF\_ITER0, SCF\_ITER1 lines. It can terminate the SCF iteration before the maximum steps (NITER0, NITER1) have been reached.

>
>TIP: Since E\_ERROR is only determined by one number (the total energy), so it can accidently reach the convergence. It might be dangerous to relied on this error to stop the SCF iteration. To avoid that, one can  reduce this value, e.g., to 1.E-8 ($eV$).

### RHO\_ERROR

|Tag|RHO\_ERROR|
| --- | --- |
|**Format**|RHO\_ERROR = value|
|**Default**|RHO\_ERROR = 0.5E-4|

The error tolerance (convergence criterion) using the SCF iteration difference between the input and output charge density. If the relative error of input and output charge density in one SCF step is less than RHO\_ERROR, the SCF iteration will be stopped.

### WG\_ERROR

|Tag|WG\_ERROR|
| --- | --- |
|**Format**|WG\_ERROR = value|
|**Default**|WG\_ERROR = 1.0E-4|

The error tolerance (convergence criterion) for the wave function conjugate gradient iterations (Hartree). This is related to SCF\_ITER0, SEC\_ITER1 lines. It can terminate the CG steps before the NLINE0, NLINE1 have been reached. This is to stop the CG iterations.

### FERMIDE

|Tag|FERMIDE|
| --- | --- |
|**Format**|FERMIDE = value|
|**Default**|FERMIDE = 0.025|

FERMIDE is the same with dE in parameter SCF\_ITER, the kT equivalent energy for Fermi-Dirac formula to calculate the electron occupations.  Default = 0.025 eV.

### MIN\_SCF\_ITER

|Tag|MIN\_SCF\_ITER|
| --- | --- |
|**Format**|MIN\_SCF\_ITER = value|
|**Default**|MIN\_SCF\_ITER = 1|

This parameter specifies the minimum number of SCF iterations steps. One can set MIN\_SCF\_ITER to a value between 2 and 8 when JOB = RELAX/MD for a more reliable result.

### MAX\_SCF\_ITER

|Tag|MAX\_SCF\_ITER|
| --- | --- |
|**Format**|MAX\_SCF\_ITER = value|
|**Default**|MAX\_SCF\_ITER = 100|

This parameter specifies the maximum number of SCF iterations steps, i.e. the NITER* in SCF\_ITER* lines. One can set MAX\_SCF\_ITER with value 1 to 1000. The SCF\_ITER* lines has higher priority than MAX\_SCF\_ITER. If the SCF\_ITER* lines are not specified, the program will use MAX\_SCF\_ITER to control the maximum number of SCF iterations steps.

### SCF\_ITER0

|Tag|SCF\_ITER0\_1/2/3/...|
| --- | --- |
|**Format**|SCF\_ITER0\_1 = NITER0\_1, NLINE0, imth, icmix, dE, Fermi-Dirac|
||SCF\_ITER0\_2 = NITER0\_2, NLINE0, imth, icmix, dE, Fermi-Dirac|
|**Default**|SCF\_ITER0\_1 = 6 4 3 0.0 0.025 1|
||SCF\_ITER0\_2 = 94 4 3 1.0 0.025 1|

These variables control the charge density self-consistent iterations for the first SCF run for JOB = SCF, RELAX, MD. For RELAX, MD, the first step SCF run (with the initial atomic positions) uses **SCF\_ITER0** lines , and subsequent steps (for moved atomic positions) uses the values of **SCF\_ITER1** lines. They are set differently because normally the first run requires much more steps. It is also used for NONSCF run, in which the ICMIX = 0 for all the NITER0 (NITER0\_1+NITER0\_2+...) lines in the following. This variable is not used for JOB = DOS.

**NITER0** is the number of SCF iterations steps. It is the sum of NITER0\_1, NITER0\_2, NITER0\_3, ... The default value for NITER0 is 100. Note the SCF iteration can be stopped before the NITER0 has been reached if the E\_ERROR has been satisfied, or the condition specified by FORCE\_RELATIVE\_ERROR has been reached. So, the stopping of SCF iteration is controlled by four parameters: **NITER0, E\_ERROR, RHO\_ERROR, FORCE\_RELATIVE\_ERROR**, whichever is satisfied first.

**NLINE0** is the number of CG line minimization steps to solve the wave functions according to $H\psi_i=\varepsilon_i\psi_i$ for a given potential (hence $H$) at each charge self-consistent step. The default value of NLINE0 is 4. Note, the CG line minimization can be stopped if the error is smaller than WG\_ERROR, or the condition specified by RHO\_RELATIVE\_ERROR is reached. So, the stopping of CG iterations is controlled by three parameters: **NLINE0, WG\_ERROR, RHO\_RELATIVE\_ERROR, whichever is satisfied first**.

These are the NITER0 lines immediately following the SCF\_ITER0 line. It describes the detail procedures for each SCF step.

**IMTH=1**, the old band-by-band CG algorithm. It should not be used unless for some special situation.

**IMTH=3**, the all band conjugate gradient method. This is the default method. We strongly recommend the use of this method.

**IMTH=2**, the DIIS method. This could be faster than IMTH = 3, but could also have stability problems. It should only be used in SCF iteration steps where the wave function is in some degree converged (e.g., not for random wave functions).

**ICMIX=0**, no charge mixing and update at this SCF step. In other word, at this step, it is a NONSCF step. For JOB = SCF, RELAX, MD, by default, for the first four SCF steps, ICMIX = 0, and ICMIX = 1 for subsequent steps. For JOB = NONSCF, for all steps, ICMIX = 0.

**ICMIX=1**, with charge mixing and update for this SCF step. Note, this is a floating point number.

Note: cne can specify something like ICMIX=1.05, as a parameter for Kerker mixing, sometime this can significantly increase the convergence speed. For most cases, ICMIX=1.00 is good enough.

**DE**: the $kT$ equivalent energy (in $eV$) for Fermi-Dirac formula to calculate the electron occupations of the eigen wave functions according to their eigen energies $\varepsilon_i$. The default value is $0.025eV$. For semiconductor, especially for defect calculation, $0.025 eV$ should be used. However, for metallic system where there are many states near the Fermi energy, one might choose a larger value, e.g., $0.1eV$ or even $0.2eV$.

**FERMI-DIRAC**: (with possible values: 0, 1, 2, 3, 4, 5). Different formulas for the Fermi-Dirac-equivalent function to calculate the wave function occupation using $\varepsilon_i$ and $dE$. These formulas are: 0, need input external files `IN.OCC` for SPIN = 1 and `IN.OCC`, `IN.OCC\_2` for SPIN = 2; 1, Fermi-Dirac; 2, Gaussian; 3,4,5 Gaussian with other prefactor polynomials. The default value is 1. However, for metallic systems, one might like to choose 2,3,4, with larger DE values.

Files IN.OCC, IN.OCC\_2 format,

```bash
1.0 1.0 1.0 0.6 0.0 0.0 0.0 ... #occupations for k-point1
1.0 1.0 1.0 0.6 0.0 0.0 0.0 ... #occupations for k-point2
```

### SCF\_ITER1

|Tag|SCF\_ITER1\_1/2/3/...|
| --- | --- |
|**Format**|SCF\_ITER1\_1 = NITER1\_1, NLINE1, imth, icmix, dE, Fermi-Dirac|
||SCF\_ITER1\_2 = NITER1\_2, NLINE1, imth, icmix, dE, Fermi-Dirac|
|**Default**|SCF\_ITER1\_1 = 40 4 3 1.0 0.025 1|

These variables control the charge density self-consistent iterations for the subsequent SCF runs for JOB = RELAX, MD. For RELAX, MD, the first step SCF run (with the initial atomic positions) uses SCF\_ITER0 lines , and subsequent steps (for moved atomic positions) uses the values of SCF\_ITER1 lines. They are set differently because normally the first run requires much more steps. It is also used for NONSCF run, in which the ICMIX = 0 for all the NITER0 (NITER0\_1+NITER0\_2+...) lines in the following. This variable is not used for JOB = DOS.

### SCF\_OUT\_FERMI\_POS

|Tag|SCF\_OUT\_FERMI\_POS|
| --- | --- |
|**Format**|SCF\_OUT\_FERMI\_POS = T / F  0/1|
|**Default**|SCF\_OUT\_FERMI\_POS = F 0|

This parameter is used to set fermi energy for calculations of insulator. If SCF\_OUT\_FERMI\_POS = T 0, output VBM  in file OUT.FERMI as the fermi energy.  If SCF\_OUT\_FERMI\_POS = T 1, output CBM  in file OUT.FERMI as the fermi energy.

### PULAY\_MIX\_OPT

|Tag|PULAY\_MIX\_OPT|
| --- | --- |
|**Format**|PULAY\_MIX\_OPT = MAX\_PULAY\_LENGTH, IFLAG\_PULAY\_WRAP, IFLAG\_OUTPUT\_PULAY, WEIGHT\_Q0\_PULAY, PENALTY\_AA\_PULAY|
|**Default**|PULAY\_MIX\_OPT = 30 1 0 1.0 0.0|

Pulay mixing method will use charge density of previous steps to precondition the input charge density of next SCF iteration. This parameter provides some options to adjust the pulay mixing.

MAX\_PULAY\_LENGTH is the maximal steps used for pulay mixing.

IFLAG\_PULAY\_WRAP indicates whether to restart pulay mixing when total pulay mxing step is great than MAX\_PULAY\_WRAP. When total pulay mixing steps is great than MAX\_PULAY\_WRAP,  if IFLAG\_PULAY\_WRAP = 0, discards all datas of previous steps and set total mixing step to zero; if IFLAG\_PULAY\_WRAP = 1, discards the datas of the one most early step, continues to do pulay mixing with the datas of previous MAX\_PULAY\_LENGTH steps. For heterostructure and low dimensional system with vacuum, MAX\_PULAY\_LENGTH = 10 and IFLAG\_PULAY\_WRAP = 0 are recommended for better SCF convergence.

If IFLAG\_OUTPUT\_PULAY = 1, output move pulay mixing infomation in REPORT. Default is 0.

WEIGHT\_Q0\_PULAY is the pulay mixing weight for charge density (in G-space) at G=0, Default is 1.0. This weight is useful for FIX\_FERMI = T. A possible setting for FIX\_FERMI = T:

PULAY\_MIX\_OPT = 100, 1, 1, 0.0001, 1.0

PENALTY\_AA\_PULAY is used to add a penalty weight to datas of previous steps, you can change the strength of the penalty by adjust this parameter. Default is 0.0, i.e. no penalty.

### PULAY\_KERK\_PARAMETERS

This includes several parameter which can be tuned to improve the SCF converged in Pulay mixing and Kerk mixing:

KERK\_AMIN = a\_min (default 0.3)

KERK\_AMIX = a\_mix (default 0.4)

KERK\_AMIX\_MAG = a\_mix\_mag (default 0.4)

KERK\_BMIX = b (default 0.5)

KERK\_BMIX\_MAG = b (default 1.e-5)

LDAU\_MIX = ldau\_mix (default 0.7)

PULAY\_WEIGHT\_SPIN = pulay\_weight\_spin (default 1.0)

PULAY\_WEIGHT\_NS = pulay\_weight\_ns (default 1.0)

In Kerk mixing, for a given $V_{in}$ and $V_{out}$ pair (comes out from the Pulay mixing), we get a new $V_{in}$ for the next SCF interation as (in G-space):

$$V_{in}(G)={0.5 a G^2 + c \over 0.5 G^2 +b} V_{out}(G)+(1-{0.5 a G^2 +c \over 0.5 G^2 + b}) V_{in}(G)$$

Here c is from the line: "SCF\_ITER0\_2 = 100,4,4,1.0x,0.025, 2", here the c=x. By default, c=0. So, one can use b to control the transition from the big G to small G region, and a is used to control the mixing rate at large G, and c is used to control the mixing rate at small G. We found that, some time using none zero c (e.g., 0.02 or 0.05) can accelerate the converges for non metallic system, and sometime even for metallic system. The default values for a and b are usually good enough.

Note, for spin=2, for the Kerk mixing, the charge densities are group into total charge and magnetic charge. The above mixing formula only applies to the total charge. For the magnetic charge, a simple mixing scheme is used, the the mixing parameter is 1 (which means the output magnetic charge density is used as the input for the next iterations). Right now, there is no input parameter to control this magnetic moment mixing.

ldau\_mix is used to control the simple  mixing parameter for the LDA+U local orbital n occupation number. Simply the $n(new)=ldau\_mix * n_{out}+ (1-ldau\_mix)* n_{in}$, here n is the local orbital occupation number matrix in LDA+U. Note, the default value for this parameter is 0.7, but sometime one can use a large ldau\_mix to accelerate the converges, e.g., even ldau\_mix=5.

pulay\_weight\_spin control the weight we used in Pulay mixing for SPIN=2 calculations. In SPIN=2, the charge densities of spin up and down are recombined into total charge density and magnetic density (the difference between up and down). The pulay\_weight\_spin control what weight we give to the magnetic density when we carry out pulay mixing. Pulay mixing is done by mixing the previous density in and out pair, trying to reduce the resulting in and out difference. When we judge how large is the in-out difference, we need to use a weight for the total charge part and spin part.

pulay\_weight\_ns is the weight factor in the pulay mixing for the case of LDA+U to place in the local orbital occupation matrix n, against the charge density.

### NONLOCAL

|Tag|NONLOCAL|
| --- | --- |
|**Format**|NONLOCAL = 1 / 2|
|**Default**|NONLOCAL = 2|

The nonlocal is the nonlocal pseudopotential implementation flag.

NONLOCAL = 1, no nonlocal potential. One has to know what he/she is doing. This usually should never be used, unless for testing purpose.

NONLOCAL = 2, the default, real space nonlocal pseudo potential implementation. It used the mask function method.

NONLOCAL = 3, g space nonlocal pseudo potential implementation. 

## Exchange-correlation tags

### XCFUNCTIONAL

|Tag|XCFUNCTIONAL|
| --- | --- |
|**Format**|XCFUNCTIONAL = LDA / PBE / HSE / PBESOL / PW91 / TPSS / SCAN / LDAWKM / LDAWKM2 / XC\_LDA\_X+XC\_LDA\_C\_PZ / XC\_GGA\_C\_PBE+XC\_GGA\_X\_PBE / XC\_HYB\_GGA\_XC\_B3LYP / XC\_MGGA\_C\_TPSS+XC\_MGGA\_X\_TPSS / SCAN + rVV10 / LDA/PBE/HSE + rVV10 / LDAWKM/LDAWKM2|
|**Default**|XCFUNCTIONAL = PBE|

XCFUNCTIONAL is used to Control the exchange-correlation functional. PWmat supports the LIBXC library for its LDA/GGA/METAGGA functional. The most commonly used functional include: LDA, PBE, HSE, PBESOL, PW91, TPSS etc. If you want to use other functional from LIBXC, you can set xcfunctional as indicated about.

If XCFUNCTIONAL = HSE, PWmat will do HSE calculation. One could use the optional parameter: HSE\_OMEGA, HSE\_ALPHA. Their default values are 0.2, 0.25 (the HSE06). Note PBE0 (the long range exchange integral) is just a special form of HSE with HES\_OMEGA parameter close to zero (see the section HSE\_OMEGA, HSE\_ALPHA). 

For hybrid-functionals in libxc except HSE and B3LYP, you should set HSE\_ALPHA and HSE\_BETA that consistent with your chosen xc functional. PWmat will set default HSE\_ALPHA and HSE\_BETA for B3LYP and HSE, but not for others.

If there are rVV10, one can optionally use RVV10\_DETAIL to specific RVV10 parameters.

The LDAWKM and LDAWKM2 options are really WKM calculations (together with JOB=SCF). In particular, what it is doing is the following: For LDAWKM:

$$H=H_{LDA} + \sum_k \lambda_k s_k (1-s_k) $$

and for LDAWKM2:

$$H=H_{LDA} + \sum_k \lambda_k s_k $$

Here $s_k$ are the occupation number of Wannier functions $\phi_k$ defined as:

$$S(k1,k2)= \sum_j <\psi_j|\phi_{k1}> <\phi_{k2}|\psi_j> occ(j) $$

Here $occ(j)$ is the occupation of the Kohn=Sham wave function $\psi_j$. Thus: $s_k=S(k,k)$.
    
The parameters $\lambda_k$ are input from the file IN.WANNIER\_PARAM. The Wannier functions in spin-up channel and spin-down channel are input from files: IN.WANNIER\_$k$.u, IN.WANNIER\_$k$.d respectively. You can check \ref{otherinput:in.wannier} for details.

For the LDAWKM, LDAWKM2 calculations, besides the usually results shown in REPORT, there are also results in OUT.WANNIER\_S (the diagonal result of s\_k) and OUT.WANNIER\_SS (the full S(k1,k2) matrix) for all the Wannier functions.

### HSE\_DETAIL

|Tag|HSE\_DETAIL|
| --- | --- |
|**Format**|HSE\_DETAIL = HSE\_MIX, MAX\_SXP, TOLHSE\_MIX, HSE\_DN, HSE\_PBE\_SCF, CHECK\_DSXH|
|**Default**|HSE\_DETAIL = 1, 1, 0.0, 6, 1, 1|

This is an optional choice when functional=HSE. It specifies all the details for a HSE SCF calculation.
   
Besides these parameters, there are other parameter which can also affect the convergence of the HSE atomic relaxation or phonon calculation, in particular the Ecutp parameter. Ecutp=4Ecut might be needed to provide sufficient accurate force. However, since increase Ecutp can significantly increase the computational time, one only use this option is Ecutp=Ecut (the default) has some problems. See the section for Ecutp.

**HSE\_MIX**, the Fock exchange kernel mixing parameter during SCF calculation (not confuse this with the HSE alpha parameter in the functional itself). It could be larger than one (e.g., 1.2). This is a bit like the charge mixing factor. Default is 1. Recommend 1 for most cases. If it is too large, it can blow up the convergence.

**MAX\_SXP**, the maximum number of Fock exchange kernel mixing terms. Numbers larger than one (e,g., 2, 3) can speed-up the HSE convergence. But each increase will cost one extra memory usage at the size of a wave function. This is like the length of Pulay mixing for charge mixing algorithm. But it is for the Fock exchange term. The default value is 1 (no Pulay mixing).

**TOLHSE\_MIX**, the tolerance for the Fock exchange term mixing for the HSE SCF iteration to stop. The default value is 0.d0. One can also use value 1.E-03. Note, in the output (screen printing), REPORT , there is one line:
```bash
update_sxp(err)(eV) = 0.2222E-02, 0.1111E-02
```

**UPDATE\_SXP** correspond to the TOLHSE\_MIX value. The first value (0.2222E-02) represents the actual Fock exchange term changes (error) after the Fock exchange term 'sxp' has been updated. The second value (0.1111-02) represent the predicted value (error) after doing Fock exchange pulay mixing when MAX\_SXP $>$ 1.

**HSE\_DN**, the number of SCF steps for each Fock exchange kernel update. Default value is 6. Recommend 3 to 10.  In our algorithm, each Fock exchange kernel update is followed by HSE\_DN steps of the SCF step (without updating the Fock exchange kernel). This HSE\_DN steps have the cost of PBE to run for each step, thus it is relatively cheap. Since each Fock exchange kernel update is expensive, this parameter is important. One wants the HSE\_DN SCF step to converge the charge density etc. following each Fock exchange kernel update. But too large HSE\_DN might not be so helpful and can still be costly to run. So, one wants to choose this parameter carefully. If one finds the HSE calculation does not converge, one might want to increase HSE\_DN. Note, the total number of SCF steps specified in SCF\_ITER0\_X, SCF\_ITER1\_X counts both the SCF update step and the Fock excahnge update steps.

**HSE\_PBE\_SCF**, the parameter define, for HSE SCF calculation, whether one wants to do one PBE calculation first. The default is 1 (do PBE SCF calculation first). If it is set to 0, that means doing the HSE calculation first, without doing a pre-PBE calculation.

**CHECK\_DSXH**, the parameter define, for HSE SCF calculation, whether check the convergency of ACE projector. The check may fail that SCF can not get the demanded accuracy(check REPORT file:ending\_scf\_reason = d\_sxp.gt.d\_sxp\_laststep ), then one can try CHECK\_DSXH=0. The default is 1.

Overall, we recommend to use HSE\_MIX=1, MAX\_SXP=1, TOLHSE\_MIX=0.0, HSE\_DN=3-6, HSE\_PBE\_SCF=1. If there are problems to converge the HSE SCF, one might want to consider MAX\_SXP=2,3, and increase HSE\_DN.

### HSE\_OMEGA

|Tag|HSE\_OMEGA|
| --- | --- |
|**Format**|HSE\_OMEGA = value|
|**Default**|HSE\_OMEGA = 0.2|

The screening parameter for HSE like hybrid functionals. Refer to J. Chem. Phys. 118, 8207 (2003) and and J. Chem. Phys. 124, 219906 (2006) for more information. PWmat support range separated hybrid functional calculations. It separate the exchange integration into long range and short range. This parameter is used to provide the range, used in the way of $\omega r$. The default value is 0.2 1/\AA. So, bigger this value, shorter the cut-off to distinguish the short range and long range.

### HSE\_ALPHA

|Tag|HSE\_ALPHA|
| --- | --- |
|**Format**|HSE\_ALPHA = value|
|**Default**|HSE\_ALPHA = 0.25|

The mixing parameter of the explicit short range Fock exchange part. The default is 0.25. Combined HSE\_OMEGA, the exchange correlation energy is:

$$E_{xc}=\alpha*E_x(FOCK,\omega)-\alpha*E_x(PBE,\omega)+E_x(PBE)+E_c$$

Here $E_x(FOCK,\omega)$ is the explicit short range Fock exchange integral with the Coulomb integration truncated by $\omega$, $E_x$(PBE) is the PBE exchange density functional energy, $E_x(PBE,\omega)$ is the short range PBE exchange density functional energy with the range defined by $\omega$.

### HSE\_BETA

|Tag|HSE\_BETA|
| --- | --- |
|**Format**|HSE\_BETA = value|
|**Default**|HSE\_BETA = 0.0|

The mixing parameter of the explicit long range Fock exchange part. The default is 0.0. PWmat provides support for a separated range hybrid calculations (all to be called HSE). The default $\beta$ is zero, thus no long range part. But one can also input a nonzero $\beta$. Sometime it is argued that, to get the correct optical properties, one should use 1/dielectric-constant as $\beta$ for bulk systems.

When $\beta$ is nonzero, the final XC functional is:

$E_{xc}=\alpha*E_x(FOCK,\omega)+\beta*(E_x(FOCK,full)-E_x(FOCK,\omega))+E_x(PBE)-\alpha*E_x(PBE,\omega)-\beta*(E_x(PBE)-E_X(PBE,\omega))+E_c$.
    
Here $E_x(FOCK,\omega)$ is the explicit short range Fock exchange integral with the Coulomb integration truncated by $\omega$, $E_x(FOCK,full)$ is the full explicit exchange integral,  $E_x$(PBE) is the PBE exchange density functional energy, $E_x(PBE,\omega)$ is the short range PBE exchange density functional energy with the range defined by $\omega$.
    
It has been argued that such a range separated HSE function can be used to replacing more advanced method like GW to provide approximated electronic structures.

### HSEMASK\_PSP

|Tag|HSEMASK\_PSP|
| --- | --- |
|**Format**|HSEMASK\_PSP1 = ampl1 size1|
||HSEMASK\_PSP2 = ampl2 size2|
|**Default**|HSEMASK\_PSP1 = 0.0 0.0|
||HSEMASK\_PSP2 = 0.0 0.0|

This is a special option for HSE, for cases where one wants to use different HSE mixing parameters for different regions (e.g., atoms) in heterostructural calculations.Since the HSE mixing parameter is very much an empirical parameter, in order to get the correct band gap at different regions (e.g., in a Si/SiO2 heterostructure), one might want to use different mixing parameters for these regions [HSEMASK](). The parameter provided here can accomplish that. 

>
>**TIP**: If you use the same pseudo for different parameter, Si for example, you can make a copy of the Si.SG15.PBE.UPF and rename it as Ge.SG15.PBE.UPF, change the element line Si to Ge in the UPF file, if you want to do MD or TDDFT calculations, plesae add an extra mass line under the element line in the pseudo UPF file, etc. mass="28.085" , if you don't do this, the software will use the mass of element Ge, that will make error result.

These are parameters used to provide an element specified HSE mixing parameter HSE\_ALPHA.The size1, size2... are in Bohr unit. The parameter can adjust the gap of HSE calculation with different  setting for different atomic types.

The explicit Fock exchange integral (including the mixing parameter) is:
    
$$ E_x(Fock)= \alpha \sum_i \sum_j o(i) o(j) \int \int \psi^*_i(r)\psi_j(r) m(r) \frac{erfc(|r-r'|\omega)}{|r-r'|} m(r') \psi_i(r')\psi^*_j(r') d^3r d^3r'$$

Here the o(i) is the orbital i occupation function, and mask function m(r) is: $m(r)=1+\sum_R ampl_R exp(-(r-R)^2/size_R^2)$, here ampl$_R$ and size$_R$ are the parameters input from HSEMASK\_PSP1,etc.

Note, the effective HSE\_ALPHA is the default HSE\_ALPHA multiplied by m(r)$^2$, which is then determined by ampl1, ampl2, etc. The ampl1 itself is not HSE\_ALPHA.

### HSE\_KPT\_TREATMENT

|Tag|HSE\_KPT\_TREATMENT|
| --- | --- |
|**Format**|HSE\_KPT\_TREATMENT = flag\_double\_grid flag\_vq\_interp vq\_interp\_range|
|**Default**|HSE\_KPT\_TREATMENT = 1 0 1.0|

This is a special option for HSE, especially for JOB=NONSCF calculations for HSE bandstructures. In HSE calculation, when bandstructure is calculated, new kpoints are provided from IN.KPT, and the wave functions from JOB=SCF calculation on a MP kpoint grid k\_prime=(NQ1,NQ2,NQ3) are used for the Fock exchange integral kernel. However, since the new kpoints from the IN.KPT are not in the k\_prime grid, some special treatments are needed to deal with this. We provided two options for this. The default option is flag\_double\_grid=1, flag\_vq\_interp=0. The default value should be good for most cases. flag\_double\_grid=1 means a special treatment for q=k-k\_prime+G. A double grid is defined on top of (within) the (NQ1,NQ2,NQ3) grid (i.e, the grid consists of lines 0, 2, 4,.. etc), the Fock exchange contribution of the k\_prime for its q on the double grid will be deleted (its corresponding Coulomb kernetl vq(q) is set to zero), while multipling the contribution from the other k\_prime points by a factor of 8/7. This can be useful, e.g., to have the correct symmetry, and make sure X point result is the same as the -X point result. Note, for JOB=SCF, one can choose flag\_double\_grid=1 or 0, although we do recommend 1 for SCF.  flag\_vq\_interp=0 means the vq interpolation is not used (vq\_interp\_range is not used in this case). 
    
For JOB=SCF, we can only use flag\_vq\_iinterp=0. In this case, a shifting of k\_prime is used. Basically, the 8 k\_prime points on the cube enclosing k are shifted to k, one at a time, so to make k points on the k\_prime grid. The 8 shifting are then linearly added using  linear interpolation weights. Note, only the kernel vq are averaged, so there is no extra FFTs. We strongly recommend this option for JOB=NONSCF for most cases. However, due to the linear interpolation, on the very small energy scale, one can have a small linear kink on the bandstructure, especially at the high symmetry point. Here, we provide another option: flag\_double\_grid=0 and flag\_vq\_interp=1. In this case, the k\_prime points are not shifted, and vq are calculated using the formula, but for q close to zero, within a cutoff qcut=dq*vq\_interp\_range (dq is the NQ1,NQ2,NQ3 grid interval size), the calculated vq (which is kind of diverging for q close to zero) is switched (using a cos function) into the Gygi's formula exxvq term. The vq\_interp\_range (usually around 1) can be used to control the range of this switching (morphing). Larger the range, smoother will be the bandstructure. Note, since we are using flag\_double\_grid=0, the energy of k point even when k=k\_prime might slightly different from the SCF value when this option is used. But when flag\_vq\_interp=0 is used, on the k\_prime point, the JOB=NONSCF bandstructure value is exactly the same as the JOB=SCF result. 

## JOB-related tags

### DOS\_DETAIL

|Tag|DOS\_DETAIL|
| --- | --- |
|**Format**|DOS\_DETAIL = IDOS\_interp, NQ1, NQ2, NQ3|
|**Default**|DOS\_DETAIL = 0|


IDOS\_interp=0 : it will not use the k-point interpolation scheme for the JOB=DOS calculation.
   
IDOS\_interp=1 : it will use the standard interpolation scheme and generate OUT.overlap\_uk.
    
IDOS\_interp=2 : it will use the second order interpolation scheme and generate OUT.overlap\_uk.2.

When IDOS\_interp=1 or 2, the interpolation borrows the method used in HSE calculation, the wave functions are first FFT into real space, then the overlap between different k-points are done. For very large systems and many k-points, it might run out of memory, one can try a smaller ECUTP / P123 to reduce the memory requirement.
    
NQ1, NQ2, NQ3 : must equal to the MP\_N123 in the last SCF or NONSCF run which generated the wave functions OUT.WG (used as IN.WG in the current DOS run).

>
>**WARNING**: When IDOS\_interp = 1 / 2, kpoint parallelization is not allowed.
{: .block-warning}

>
>**WARNING**: When IDOS\_interp = 2, NQ1,NQ2, NQ3 needs to be greater than or equal to 4.
{: .block-warning}

### NUM\_DOS\_GRID

|Tag|NUM\_DOS\_GRID|
| --- | --- |
|**Format**|NUM\_DOS\_GRID = num|
|**Default**|NUM\_DOS\_GRID = 4000|

This is the number of energy grid points when calculating DOS, the default is 1500. This grid is within the window [$E_{min}$, $E_{max}$], and the $E_{min}$, $E_{max}$ are determined by the minimum and maximum eigen energies.

### DOS\_GAUSSIAN\_BROADENING

|Tag|DOS\_GAUSSIAN\_BROADENING|
| --- | --- |
|**Format**|DOS\_GAUSSIAN\_BROADENING = value|
|**Default**|DOS\_GAUSSIAN\_BROADENING = 0.05|

This is the gaussian broadening of DOS plot for JOB=DOS, unit is eV.

### RELAX\_DETAIL

|Tag|RELAX\_DETAIL|
| --- | --- |
|**Format**|RELAX\_DETAIL = IMTH, NSTEP, FORCE\_TOL, ISTRESS, TOL\_STRESS, TOL\_LINECORRECTION|
|**Default**|RELAX\_DETAIL = 1, 200, 0.02, 0, 0, -0.001 (ACCURACY = NORM)|
||RELAX\_DETAIL = 1, 200, 0.01, 0, 0, -0.001 (ACCURACY = HIGH / VERYHIGH)|
||RELAX\_DETAIL = 1, 200, 0.03, 0, 0, -0.001 (XCFUNCTIONAL = HSE)|

This is an optional line for ``JOB = RELAX''. It controls the atomic relaxation steps. Note PWmat1.5+ can relax the lattice vectors.

**IMTH**, the method of relaxation. The default is 1, conjugated gradient. Other options include 2, BFGS method, 3, steepest decent (this is mostly for JOB=NEB), 4, Preconditioned Conjugate Gradient (PCG), experimental feature, See \ref{tag:vffdetail} before you use this, 5, Limited-memory BFGS method, 6, FIRE: Fast Inertial Relaxation Engine.

Please note that: all imths (1,2,3,4,5,6) can be used for atomic relaxation,  but only imth=1,5,6 can be used for cell relaxation. In addition, some options can be set for the optimizers in file IN.RELAXOPT(need to set IN.RELAXOPT=T in etot.input), see section \ref{otherinput:in.relaxopt} for more details.

**NSTEP**, the maximum number of relaxation steps. Each total energy calculation is one step, i.e., it counts the steps inside the line minimization in the total steps. In another word, NSTEP is more (at least twice) than the CG steps.

**FORCE\_TOL**, the force tolerance for the maximal residual force. If the maximum force is less than FORCE\_TOL, the relaxation will stop.

**ISTRESS**, controls whether to relax the lattice vectors. If ISTRESS=0 (or the last two number do not exist), the lattice will not be relaxed. If ISTRESS=1, PWmat will relax lattice vectors. One can add external stress tensor by setting STRESS\_EXTERNAL or PTENSOR\_EXTERNAL in file atom.config, and external pressure by setting PSTRESS\_EXTERNAL in file IN.RELAXOPT, the latter need to set IN.RELAXOPT=T in etot.input.
    
If you have set STRESS\_EXTERNAL or PTENSOR\_EXTERNAL, make sure the settings are consistent with sysmetry operations you have used (IN.SYMM=T) or generated by MP\_N123, if not you should turn off the symmetry operations. Check MP\_N123 for details about symmetry.

If ISTRESS=1 and XCFUNCTIONAL=HSE, one should set RELAX\_HSE as follows:

```bash
    RELAX\_HSE = 0   0.0 2
```
    
Basic settings for cell relax with xcfunctional=pbe :

```bash
    JOB = RELAX
    XCFUNCTIONAL = PBE
    RELAX_DETAIL = 1 1000 0.02 1 0.05
    ECUT = 70 # maybe 1.4 * Ecut_default
    ECUT2 = 280
```

Basic settings for cell relax with xcfunctional=hse

```bash
    JOB = RELAX
    XCFUNCTIONAL = HSE
    RELAX_DETAIL = 1 1000 0.03 1 0.05
    ECUT = 70 # maybe 1.4 * Ecut_default
    ECUT2 = 280
    RELAX_HSE    =    0   0.50000E-01     2
```

**TOL\_STRESS**, the stress tolerance for the maximal residual stress. (here it is defined as $\partial{Etot}/\partial{STRAIN}/Natom$, Etot is the energy of the whole system (not the energy of unit volume)).

**TOL\_LINECORRECTION**, the energy tolerance for the line minimization alone one search direction. When we see Etot(step) becomes a linear line, perhaps we should use a smaller value to have more correction steps.  However, for more accurate relaxations, we can set a smaller value for this parameter, so the relaxation can continue. We can turn off the energy tolerance checking for the line minimization by setting TOL\_LINECORRECTION < 0. Some time energy is not that much accurate, the energy tolerance checking is not reliable. And the default value is -0.001.

The JOB = RELAX will output a RELAXSTEPS and MOVEMENT files. While RELAXSTEPS gives a summary of the steps, MOVEMENT records the atomic positions and lattice vectors for all the steps.

**SPECIAL FORCE ON EACH ATOM**:
In the RELAX calculation, one can add atom specified external force on each atom. This is done by using IN.EXT\_FORCE=T in etot.input. In that case, A file called IN.EXT\_FORCE needs to be provided, it will give the external force on each atom during the atomic relaxation. Please see the section IN.EXT\_FORCE for more details.

**SOME DISCUSSION ABOUT RELAXATION**:

Atomic relaxation is one of the most used feature in DFT calculations. One has to balance the speed with the stability. Here, we have implemented 6 different methods. For most common problem, we suggest to use imth=1 (conjugate gradient,CG).

For very large system, to accelerate the convergence, one can test the use of imth=4, which is the VFF accelerated relaxation. In the best case (e.g., very large systems), imth=4 can speed up imth=1 by a factor of 10. However, one might want to test its stability. One can also used imth=4 (the BFGS method), sometime it is faster than the CG method. Another new method is imth=6, it uses a molecular dynamics but with a friction term, so the system will eventually relax to a local minimum. This can be stable, but one needs to test the parameter FIRE\_DT in IN.RELAXOPT. 

Lastely, if all these methods are unstable,  one can always use imth=3, and use a very small maximum step by setting RELAX\_MAXMOVE in IN.RELAXOPT. This might take a long time, but it should be stable if sufficient small RELAX\_MAXMOVE is used.

When doing atomic relaxation, one must be mindful of the egghead problem, which is the artificial forces caused by the real space lattice. One can remove this problem by setting Ecut2=4Ecut, and Ecut2L=4Ecut2. Most likely, the second condition Ecut2L=4Ecut2 might not be necessary, or one can use JOB=EGGFIT and EGG\_DETAIL to remove the egghead effect without using Ecut2L=4Ecut2. Unfortunately, one might always need to use Ecut2=4Ecut.

To check the relaxation convergence, one should always check the energy as a function of iteration steps in RELAXSTEPS. Note, that, some steps might be trial steps, so they are not so important (you might see some spike). The important one is the "NEW" step energy.

If funcitonal=HSE is used for atomic relaxation, some special algorithms are used for its acceleration. Please check RELAX\_HSE for details.

Finally, if cell lattices are relaxed, great care must be taken. During the lattice relaxation, the number of plane wave basis (the G-vectors within the Ecut sphere) is not changed. So, at the end of the relaxation, the plane wave set (it  might no longer be a sphere) might not correspond to the Ecut sphere. So, if one uses the Ecut to run it again, a new plane wave set based on the Ecut sphere will be chosen, and the energy and the stress might be changed. So, either one needs to do this multiple times, or one needs to use a rather large Ecut, so the change of basis will have minimum effect. One can use the STRESS\_CORR to mitigate (compensate) this problem in some degree, but cannot really remove it completely. This will be particularly problematic if the volume change is very large. In some cases, it might be more reliable to do the cell relaxation by hand, specially if only one degree of freedom is used. Note, one can use stress\_mask in the atom.config to choose what lattice components can be changed. One can also used imov in atom.config to determine which atom, and which x,y,z direction can be moved for internal atomic relaxation.

### RELAX\_HSE

|Tag|RELAX\_HSE|
| --- | --- |
|**Format**|RELAX\_HSE = NUM\_LDA, FACT\_HSELDA, LDA\_PBE|
|**Default**|RELAX\_HSE = 20, 0.05, 2|

This is an optional line for ``JOB = RELAX'' when XCFUNTIONAL=HSE. It uses special techniques to accelerate the atomic relaxation under HSE. Currently, this option only works for conjugated gradient atomic relaxation as defined in {RELAX\_DETAILS}. In this option, the additional LDA or PBE atomic relaxations are used as preconditioner for the HSE relaxation.

**NUM\_LDA** is the maximum number of relaxation steps for the LDA/PBE preconditioner run. It uses the LDA/PBE atomic relaxation as a preconditioner for the HSE atomic relaxation.  If {NUM\_LDA=0}, then no precondition is used, it is the plain CG atomic relaxation based on HSE. The default is 20.

**FACT\_HSELDA** is the prefactor to stop the LDA/PBE relaxation: if LDA/PBE force is less than FACT\_HSELDA multiplied the HSE force, then stop. The default is 0.05.

**LDA\_PBE** is the indicator for LDA or PBE functional used for the atomic relaxation to find the preconditioner of HSE relaxation. If LDA\_PBE=1, use LDA; LDA\_PBE=2, use PBE. One should use the xcfunctional which is closest to HSE. The default value is 2.

Comments: before one use this option (NUM\_LDA $>$0), one better make sure the PBE, or LDA atomic relaxation of the system is smooth. So, one might want to use Ecut2=4*Ecut.

### VFF\_DETAIL

|Tag|VFF\_DETAIL|
| --- | --- |
|**Format**|VFF\_DETAIL = FF\_IMTH, FF\_NSTEP, FF\_FORCE\_TOL, K\_BOND, K\_ANGLE, K\_DIHEDRAL, K\_SHIFT|
|**Default**|VFF\_DETAIL = 1, 500, 0.01, 30.0, 5.0, 0.0, 1.5|

Note, if you want to use PCG(IMTH=4) method to accelerate the relaxation, you should know some basic concepts about force field theory. The PCG method supports any systems such as molecular, metallic, semiconductor and it has a well optimized energy function to the metallic systems like Al, Ni, Au, Cu, Ag, Pt, Ir, Pd, Rh, La, Ce, Mg, Ca, Sr. This method has a high stability for molecular systems, semiconductors and gives high speedup factors for molecule absorbed on surface cases.

**FF\_IMTH**: the optimization algorithm used in force field relaxation.

**FF\_NSTEP**: the maximum optimized steps in force field relaxation.

**FF\_FORCE\_TOL**: the tolerance of force in force field relaxation.

**K\_BOND, K\_ANGLE, K\_DIHEDRAL, K\_SHIFT**: the force constant for bond, angle, dihedral, shift terms.
$$
        E_{tot} = k_b(b-b_0)^2 + k_a(\theta-\theta_0)^2 + k_d(\phi-\phi_0)^2 + k_s(r-r_0)^2
        $$

The defaults is \textbf{VFF\_DETAIL = 1, 500, 0.01, 30.0, 5.0, 0.0, 1.5}. Note, it is not recommended to modify these parameters unless you are an expert on force field. If there is a metallic area in your system, that is, some metallic atoms have only metallic atoms neighbors and only metallic bond, you have to set additional parameters in the `atom.config' file like this:
```bash
29    0.22    0.11    0.06   1  1  1  1  1
 1    0.44    0.49    0.39   1  1  1  1  0
```
Here the Copper atom(Z=29) has five integers after coordinates, first three integers determine whether the atom will move along particular direction, the forth integer number has some meaning for DOS calculation and the last integer number (the ninth column) determine whether the atom is in the metallic area. As we know, the hydrogen atom is not a metallic atom, we set the last integer as 0, while Copper is metallic, so we have set it to 1. Note, you can set Copper to 0, so it will be dealt as a covalent bond element. Usually we only set it to 1 when there is a large piece of metal.

### MD\_DETAIL

|Tag|MD\_DETAIL|
| --- | --- |
|**Format**|MD\_DETAIL = MD, MSTEP, DT, TEMP1, TEMP2|
|**Default**|None|

This line is required when JOB=MD, or JOB=TDDFT, or JOB=NAMD. There is no default values, hence must be input by hand.

**MD** : the method of MD algorithm

|MD|Method|
| --- | --- |
|1|Verlet (NVE)|
|2|Nose-Hoover (NVT)|
|3|Langevin (NVT)|
|4|Constant pressure Langevin dynamics (NPT)|
|5|Constant pressure Nose-Hoover dynamics (NPT)|
|6|Berendsen dynamics (NVT)|
|7|Constant pressure Berendsen dynamics (NPT)|
|8|Multi-Scale Shock Technique (MSST)|
|11|Verlet (NVE), continue run|
|22|Nose-Hoover (NVT), continue run|
|33|Langevin (NVT), continue run|
|44|Constant pressure Langevin dynamics (NPT), continue run|
|55|Constant pressure Nose-Hoover dynamics (NPT), continue run|
|66|Berendsen dynamics (NVT), continue run|
|77|Constant pressure Berendsen dynamics (NPT), continue run|
|88|Multi-Scale Shock Technique (MSST), continue run|
|100|Calculating multiple configurations stored in IN.MOVEMENT|
|101|Calculating multiple configurations stored in IN.MOVEMENT|


Verlet is for NVE (fixed number of atom N, fixed volume V, and fixed total energy E), Langevin and Nose-Hoover are for NVT (fixed number of atom N, fixed volume, and fixed temperature T), and Constant pressure Langevin or Nose-Hoover dynamics are for NPT (fixed number of atom N, fixed pressure P, and fixed temperature T). Currently, we do not have NPE. We also provide MSST (Multi-Scale Shock Technique) to simulate a compressive shock wave passing over the system.

One can also set: MD=11,22,33,44,55,66,77,88 which means the continue run of MD following the previous runs.

**MSTEP** : the number of MD steps.

**DT** : the time length for each MD step (in the unit of $fs$, $1fs =1\times10^{-15}s$). Note, usually, with H atoms, dt should be $1fs$, and with heavier atoms, dt could be $2fs$. However, for rt-TDDFT run, dt should be much smaller, like 0.1 fs to 0.2 fs.

**TEMP1** and **TEMP2** : the beginning and final temperature (in $Kelvin$). When there is no velocity session in the atom.config file, the TEMP1 will be used to randomly generated an initial velocity (the initial kinetic energy is generated as twice the 0.5*K*T, with the expectation that half of its energy will be converted into potential energy). So, in the simulation, the istep=1 temperature will be 2*TEMP1.

TEMP2 will not be used for MD=1, or 11 (NVE) (but still, it should be there as a place holder). During the MD, the program will adjust the temperature linearly, let it goes from TEMP1 to TEMP2. For MD=3,5 (Langevin), one can use a LANGEVIN\_ATOMFACT\_TG section in the atom.config file to specify a local atomic specified temperature (and Langevin parameter gamma). In that case, the desired atomic temperature at a given time equals the global desired temperature calculated from TEMP1 to TEMP2, then multiplied by the atomic scaling factor specified in the LANGEVIN\_ATOMFACT\_TG.

For method 1-8, one can use file IN.MDOPT to set detailed parameters by setting IN.MDOPT=T, all the parameters will be written in file OUT.MDOPT. Please refer to section \ref{otherinput:in.mdopt} for details.

In the Berendsen method (MD=6,7,66,77), the kinetic energy is scaled at every MD step as (1+(Tdesired/Tcurrent-1)*dt/tau). For MD=7,77 (NPT), the cell box is scaled at every MD step as:

$$
cell(i1,i2) = cell(i1,i2) * ( 1+(press(i1,i2) - press\_ext(i1,i2)) * stress\_mask(i1,i2) * dt / tauP ).
$$

Here press, press\_ext are pressures in the unit of $eV/Angstrom^3$. Note, press(3,3) equals the stress(3,3)/volume, so it is a 3x3 tensor. The external pressure $press\_ext(i,j)=delta\_{i,j} MD\_NPT\_PEXT\_XYZ(i)$ . Note, this is input from IN.MDOPT, not the stress from the atom.config. If MD\_NPT\_PEXT\_XYZ is not specified, then MD\_NPT\_PEXT\_XYZ(:) = MD\_NPT\_PEXT. If both MD\_NPT\_PEXT\_XYZ and MD\_NPT\_PEXT are not speficied, then the external pressure is zero. The stress\_mask is from atom.config. The detault value is stress\_mask(i1,i2)=1 (for every element of the matrix).

One can always set MD\_SEED and MD\_AVET\_TIMEINTERVAL. For MD=4 or 5, one can check the internal pressure in file MDSTEPS and MOVEMENT.

In the MSST method(MD=8,88), you always need to set MD\_TAU\_CELL, MD\_MSST\_VS, MD\_MSST\_DIR. MD\_TAU\_CELL will set the masslike parameter for the simulation cell size, i.e. Q in paper \cite{Reed}, but MD\_TAU\_CELL itself is the time to arrive equilibrium. MD\_MSST\_VS is the shock speed $v_s$ in paper \cite{Reed}. MD\_MSST\_DIR is the shock direction, its value can be 0, 1 or 2 for in x, y or z direction. IN file MOVEMENTS you can check MD\_MSST\_INFO for additional outputs. IN file MDSTEPS there adds a new column "V/V0", which shows the change of vulome. In the process of MSST, if MD\_MSST\_VS is kind of large in some way and the size of box changes alot, PWmat will hard to converge or evan crash. You need to reduce the time step, i.e. use smaller DT in MD\_DETAIL.

When MD=11/22/33/44/55/66/77/88, it is a continue run for Verlet/Langevin/Nose-Hoover constant pressure, Langevin/constant pressure Nose-HooverNPT, Berendsen NVT, Berendsen NPT respectively. In these cases, the atom.config file should include the velocity section. Note, for JOB=MD, if there is velocity section in atom.config, the velocity will be used, and there is no initial scaling of the velocity using temperature TEMP1.

**SPECIAL LV**:  There is a special feature for LV dynamics (either 3 or 4). In this special feature, we can specify the desired temperature for each atom. We can also specify the desired GAMMA value (the MD\_LV\_GAMMA) for each atom. In another words, you can make one atom very hot, and let the temperature decaying from this atom. The temperatue decaying length will be controlled by the MD\_LV\_GAMMA. Smaller this value, decaying length will be longer, i.e., more graduate. The atom specific temperature and Gamma are controlled by scaling factors, they are specified in the atom.config file, with a special session called: "LANGEVIN\_ATOMFACT\_TG". They have the following formats (see section \ref{inputfile:atomconfig}):

```bash
LANGEVIN_ATOMFACT_TG
30  0.5   1.0
30  0.2   1.2
..............
atom scaleT  scaleG

```

Here the desired temperature for one atom equals the original desired temperature specified by temp1, temp2 and the steps, then multiplied by scaleT(iatom). The Gamma for one atom equals the MD\_LV\_GAMMA specified in the above table, or its default value, multipled by scaleG(iatom).

**SPECIAL FORCE**: In all the MD calculation, one can add atom specified external force on each atom. This is done by using IN.EXT\_FORCE=T in etot.input. In that case, A file called IN.EXT\_FORCE needs to be provided, it will give the external force on each atom during the molecular dynamics, or during atomic relaxation. Please see the section IN.EXT\_FORCE for more details.


**MD=100,101**: for these choices, instead of doing an actual MD simulation following the atomic forces, the multiple configurations stored in a file IN.MOVEMENT will be calculated one after another, so the forces will not really be used, and the trajectory follows the one in IN.MOVEMENT. Note, in the running directory, a IN.MOVEMENT file need to be provided. This is mostly used for some special purposes, for example for machine learning force field development, while the IN.MOVEMENT is generated by force field. For MD=101, the charge and wave function interpolation is turned off. This might be useful if the configures in IN.MOVEMENT change dramatically from one frame to another frame. While the format in IN.MOVEMENT is the same as in the output MOVEMENT, it must provide a header.

The formate of IN.MOVEMENT, 

```bash
nstep,nskip1,nskip2,nkip3,njump
64, atoms,Iteration .....
...
Lattice vector (Angstrom)
0.1130000000E+02    0.0000000000E+00    0.0000000000E+00
0.0000000000E+00    0.1130000000E+02    0.0000000000E+00
0.0000000000E+00    0.0000000000E+00    0.1130000000E+02
Position (normalized), move_x, move_y, move_z
31         0.99973    0.99973    0.99973     1  1  1
31         0.99985    0.24985    0.24985     1  1  1
.......
.......
```

The code used to read the IN.MOVEMENT file is as following,

```fortran
do istep=1,nstep
do ii=1,njump
do i=1,nskip1
read(IN.MOVEMENT,*)
enddo
read(IN.MOVEMENT,*) AL(1,1),AL(2,1),AL(3,1)
read(IN.MOVEMENT,*) AL(1,1),AL(2,1),AL(3,1)
read(IN.MOVEMENT,*) AL(1,1),AL(2,1),AL(3,1)
do i=1,nskip2
read(IN.MOVEMENT,*)
enddo
do i=1,natom
read(IN.MOVEMENT,*) iat(i),x1(i),x2(i),x3(i)
enddo
do i=1,nskip3
read(IN.MOVEMENT,*)
enddo
enddo ! ii=1,njump
enddo ! istep=1,nstep
```

The nskipt1,nskipt2,nskip3 are the skips of lines in different segment of the IN.MOVEMENT file.
The njump indicates whether you like to calculate every configuration (njump=1), or you like to
jump over some configurations, e.g., njump>1. nstep is the total number of steps to calculate.
You might need to check  the IN.MOVEMENT file to determine the nskip1,nskip2,nskip3. nskip2 is
usually 1.

### MD\_VV\_SCALE

|Tag|MD\_VV\_SCALE|
| --- | --- |
|**Format**|MD\_VV\_SCALE = NSTEP|
|**Default**|MD\_VV\_SCALE = 100|

To scale the kinetic energy in Verlet MD (for JOB=MD, iMD=1/11) for every NSTEP steps, so the total energy is conserved. The default NSTEP is 100. This is used for enforce the total energy conservation for Verlet. Note, in TDDFT, NAMD, or some MD when there is external potential or electric field, the total energy is not supposed to be conserved. In those case, please set MD\_VV\_SCALE to a very large number, so it will never be used. The default value is MD\_VV\_SCALE=100 for MD, and not used for TDDFT and NAMD.

### TDDFT\_DETAIL

|Tag|TDDFT\_DETAIL|
| --- | --- |
|**Format**|TDDFT\_DETAIL = $m_1$, $m_2$, mstate|
|**Default**|TDDFT\_DETAIL = 1, NUM\_BAND, NUM\_BAND|

This will be read in when JOB = TDDFT. Note if mstate=-1, this is for TDDFT\_NOB calculation, see below.

Note, when JOB=TDDFT, besides TDDFT\_DETAIL, it also reads in parameters from MD\_DETAIL, and optionally from TDDFT\_SPACE, TDDFT\_TIME, TDDFT\_STIME.

In the TDDFT calculation, we expand the time dependent electron orbital $\psi_j(t)$ in terms of the adiabatic eigenstates $\phi_i(t)$.

$$
\psi_j(t)=\sum\limits_{i}C_{ji}\phi_i(t)
$$

There will be $mstate$ electron wavefunctions [j=1,mstate] which will be occupied by their occupation number o(j) and described by $\psi_j(t)$.

However, for the first $m_1-1$ orbital ($j=1,m_1-1$), $\psi_j(t)$ is just $\phi_i(t)$, i.e., these $m_1$ states are fully occupied like in Born-Oppenheimer MD, no electron excitation:

$$
\psi_j(t)=\phi_j(t), j=1,m_1-1
$$

For the next $j=m_1, mstate$ state (so, mstate include the [1,$m_1-1$] state!), we will expand the $\psi_j(t)$ in the $\phi_i(t)$ window of $i=m_1,m_2$:

$$
\psi_j(t)=\sum_i C_{ji}(t)\phi_i(t), j=m1,mstate;i=m1,m2
$$

Thus, in total, the expansion window is [$m_1,m_2$], and the total number of time dependent orbital is: mstate (they will be occupied by o(j), so mstate can be larger than the NUM\_ELECTRON/IPSIN. However, within the mstate, the first $m_1-1$ states are fully occupied, and just equal to the adiabatic eigen states, the next $mstate-m_1$ state will be expanded using adiabatic state within the window of [$m_1,m_2$], and their occupation might follow the Fermi-Dirac rule, or to be input by IN.OCC/IN.OCC\_2 (see Appendix \ref{appendix:tddft}). The initial $C_{ji}$ can also be input from IN.CC/IN.CC\_2 (see Appendix \ref{appendix:tddft}).

|||
| --- | --- |
|**[m1,m2]**|Adiabatic window $\phi_{i,i=m1,m2}$. The $[1,m1-1]$ will always be occupied by the first $\psi_{j,j=1,m1-1}$ states. $m2 \in [m1,NUM\_BAND]$ , usually $m2$ is smaller than $NUM\_BAND$ by a few states, because the last few states maybe not converge well.|
|**[1,mstate]**|Wavefunction index. $\psi_{j,j=1,mstate}$. $mstate\in [m1,m2]$|

The choice of m2 is important for the physical correctness of the TDDFT simulations. The choice might depend on the physical problems at hand. Larger the m2, more accurate will be the simulation, but it can also cost more time to calculate. Typically, from mstate to m2, one should include all the possible electron excitations. For example, for a light excitation, if the hot electron can be 1-2 eV above the bottom of conduction band, then m2 should be choosen to include all these bands. Sometime m2 can be twice as mstate. But tests are needed to determine this. All depend on how high the electron can be excited to the conduction band.

If mstate=-1, this is a special case, for TDDFT\_NOB calculation. NOB stands for natural orbital branching. In this case, we need another line, which is:

```bash
TDDFT_NOB = iseed, S_c, tau, temp2, kin_scale, select_opt
```

iseed (negative integer) is a random number seed for the stochastic NOB calculation. S\_c is the cut-off entropy for the branching. tau is the dephasing time (in fs). If tau is negative, then IN.BOLTZMANN\_TAU will be used to input tau(i) for each state, and the tau$_{ij}$ will be determined from $\sqrt{tau(i)tau(j)}$. temp2 is the  temperature for the case of kin\_scale=2. kin\_scale=1,2,3 are three different ways for kinetic energy scaling after branching, similar to the flag\_scale=1,2,3 in the TDDFT\_BOLTZMANN flag.

kin\_scale=1 will scale all the atom's velocity uniformly to conserve the total energy.

kin\_scale=2 will scale all the atom's velocity uniformly to keep the temperature at temp2. As a result this method will not conserve the total energy, and the total energy will usually graduately decrease.

kin\_scale=3 will be the standard way to scale the kinetic energy in the transition degree of freedom, and conserve the total energy.

Select\_opt=1,2,3 is the option for branching algorithm.

1: for using the transition eigen energy difference and the current temperature (derived from the current kinetic energy) with an Boltzmann factor to determine the natural orbital branching probability to restore the detailed balance.
    
2: for using the SCF total energy difference after a trial branching and the current temperature to determine the natural orbital branching probability.

3: this will also use the actual SCF total energy, instead of eigen energies, to determine whether a branching is allowed. However, instead of using a temperature and an Boltzamn factor, in this scheme, the transition degree of freedom is first determines, and whether this degree of freedom has enough kinetic energy to compensate the SCF total energy increase, is used to determine whether one particular branching is allowed.  Note, for select\_opt=3, one must also choose kin\_scale=3.

Note, kin\_scale=3, select\_opt=3 will be the standard way of doing statistical branching.

### TDDFT\_SPACE

|Tag|TDDFT\_SPACE|
| --- | --- |
|**Format**|TDDFT\_SPACE = itype\_space, N, a(1), ..., a(N)|
|**Default**|TDDFT\_SPACE = 0|

This controls the real space Vext\_tddft(r). Vext\_tddft(r) refers to the  external potential in real space for tddft calculation.

|itype\_space|Description|
| --- | --- |
|0|No external input term.|
|1|Read vext\_tddft from file IN.VEXT\_TDDFT(all capital, same format as in IN.VEXT).|
|2|$Vext\_tddft(r)=(x-a(1))a(4)+(x-a(1))^2a(5)+(y-a(2))a(6)+(y-a(2))^2a(7)+(z-a(3))a(8)+(z-a(3))^2a(9)$, a(1),a(2),a(3) in fractional coordinates, a(4)-a(8) in unit of Hartree/Bohr. output file OUT.VEXT\_TDDFT.|
|3|$Vext\_tddft(r)=a(4)e^{-[(x-a(1))^2+(y-a(2))^2+(z-a(3))^2]/a(5)^2}$, a(1),a(2),a(3) in fractional coordinates, a(4) in unit of Hartree, a(5) in unit of Bohr. output file OUT.VEXT\_TDDFT.|
|-1|Not use real space format, but use G-space,it wil use IN.A\_FIELD|

The `IN.VEXT\_TDDFT' file can be copied from other TDDFT calculation output file 'OUT.VEXT\_TDDFT', or generated by utility programs \hyperref[util:addfield]{add\_field.x}.

### TDDFT\_TIME

|Tag|TDDFT\_TIME|
| --- | --- |
|**Format**|TDDFT\_TIME = itype\_time, N, b(1), ..., b(N)|
|**Default**|TDDFT\_TIME = 0|

This is used to control the time dimension of the external function fTDDFT(i).

|itype\_time|Description|
| --- | --- |
|0|$ftddft(t)=1.0$|
|1|read in $ftddft(i)$ from IN.TDDFT\_TIME|
|2|$ftddft(t)=b(1)e^{-(t-b(2))^2/b(3)^2)}\sin(b(4)t+b(5))$.  $b(2)$,$b(3)$ in unit of $fs$; $b(4)$ in unit of rad/fs unit, $b(5)$ in unit of rad; $b(1)$ no unit. output file OUT.TDDFT\_TIME|
|22|$ftddft(t)=\int^t_0 [b(1)e^{-(t-b(2))^2/b(3)^2)}\sin(b(4)t+b(5))] dt$.  $b(2)$,$b(3)$ in unit of $fs$; $b(4)$ in unit of $rad/fs$, $b(5)$ in unit of rad; $b(1)$ no unit. output file OUT.TDDFT\_TIME|

File IN.TDDFT\_TIME format,

```bash
0 ftddft(0)
1 ftddft(1)
...
N ftddft(N)

```

For TDDFT Hamiltonian, we have,

|Option|Description|
| --- | --- |
|$\ne -1$ |$H(t)=H_0+Vext\_tddft(r)ftddft(t)$|
|-1|$H(t)=-1/2(\nabla_x+i A_x*ftddft(t))^2-1/2(\nabla_y+i A_y*ftddft(t))^2-1/2(\nabla_z+i A_z*ftddft(t))^2$|

### TDDFT\_BOLTZMANN

|Tag|TDDFT\_BOLTZMANN|
| --- | --- |
|**Format**|TDDFT\_BOLTZMANN = flag\_b, flag\_scale, temp, tau, istep\_start(opt), nstep\_CG(opt)|
|**Default**|TDDFT\_BOLTZMANN = 0|

This line controls whether to introduce Boltzmann factor in order to keep the correct detailed balance between two adiabatic states $\phi_{i1}(t)$ and $\phi_{i2}(t)$. This goes beyond the usual Ehrenfest dynamics. In the Ehrenfest dynamics, the electronic system will be over heated due to the lack of detailed balance (which means the transition from the lower energy adiabatic state $\phi_{i1}(t)$ (with eigen energy $E_{i1}(t)$) to higher energy adiabatic state $\phi_{i2}(t)$ (with eigen energy $E_{i2}(t)$) is a factor of $exp(-(E_{i2}-E_{i1})/kT)$ smaller than the transition from $\phi_{i2}(t)$ to $\phi_{i1}(t)$. This suppression of the up-lifting transition can be realized by adding this Boltzmann factor. Adding this Boltzmann factor is critical in order to have the proper hot electron cooling. However, an dephasing time tau has to be used with an special algorithm when adding this Boltzmann factor, otherwise the cooling will be too fast (e.g., an energy conservation between the transition states and an phonon mode in the weak coupling regime will not be satisfied). A special algorithm is implemented in PWmat to properly take into account the tau. 

In a way, this is like the Tully's surface hopping, but without the stochastic feature in the dynamics. Compared with Tully's algorithm, it has more correct dephasing behavior. Compared with wave function collapsing behavior, it can have more proper treatment for tau. However, this is a meanfield treatment for the nuclear trajectory. Unlike the stochastic potential energy surface hopping, or wave function collapsing treatment, the current meanfield treatment does not provide a branching for the trajectory (thus might not be good if you like to calculate the probability for  different chemical reaction, etc). The Boltzmann factor is not applied to each individual electron state $\psi_j(t)$ (and its collapsing), instead, it is applied to the occupation of adiabatic state $\phi_i(t)$ (collectively for all $\{ \psi_j(t) \}$),  as a result, it has a property of unitary rotational invariance among the group of $\psi_j(t)$.

**flag\_b** : indicates whether to use the Boltzmann factor: 0 means no Boltzmann factor (Ehrenfest dynamics); 1 means with Boltzmann factor. Note, one has to detemine according to the physics, whether to use this Boltzmann factor. For example, if the dynamics does not involve phonon, and it is very short (e.g., using rt-TDDFT to simulate the light absorption), then one might not want to include the Boltzmann factor. One thing needs to be mindful: for long time rt-TDDFT simulations, the use of Boltzmann factor can make the simulation more stable, since all the states are moving towards equilibrium.

**flag\_scale** : indicates what method to be used to scale the kinetic energy. When using the Boltzmann factor, the total energy will not be conserved without the rescaling of the velocity of the nuclei (usually, the total energy will decrease due to the electron cooling). There are different ways to deal with this. flag\_scale=0, this means without scaling the kinetic energy, as a result, the total energy will decrease with time. flag\_scale=1, this means the velocities of all the nuclei will be rescaled with a uniform scaling factor, so the total energy will be an constant (thus the lost electronic energy will be given to the phonon movement). For example, this will be useful for isolated molecule. In this case, if there is an initial hot electron, the temperature of the system will likely increase with time. Note, the Boltzmann factor $exp(-(E_{i2}-E_{i1})/kT)$ dynamically depends on this temperature T(t). flag\_scale=2,this means the velocities of all the nuclear will be scaled to keep the temperature (the total kinetic energy) to a constant, specified by $temp$. In this case, the temperature in the $MD\_DETAIL$ will not be used, but the $temp$ specified here will be used. This might be a good approximation if the studied system is embedded in a thermal bath, which is always kept in a constant temperature. flag\_scale=3, this is like $flag\_scale=1$, where the kinetic energy is modified to keep the total energy conserved. However, instead of uniformly scale the velocity of all the atoms by a constant factor, here the electron-phonon coupling constant is used to specify a force of a given atom, which is then used to be added to the current velocity to conserve the total energy. This is more rigorous treatment for what phonon degree of freedom to give the extra kinetic energy for. More specifically, let's use TCD(i1,i2) to indicate the Boltzmann factor introduced  modification of the density matrix, which represents the charge transfer between adiabatic states i1 and i2 (this modification maintains the detailed balance between adiabatic states i1 and i2, but also cause the violation of energy conservation).

We also define $D(i1,i2)=\sum_j C*(i1,j) C(i2,j) o(j)$, here C(i,j) is the expansion coefficient of electron wave function $\psi_j(t)$ on adiabatic state $\phi_i(t)$, i.e, $\psi_j(t)=\sum_i C(i,j) \phi_i(t)$, and o(j) is the occupation (not change with time) of electron state $\psi_j(t)$. Then the extra atomic force due to the TCD(i1,i2) between adiabatic states $\phi_{i1}(t)$ and $\phi_{i2}(t)$ will be: $F(R)= \sum_{i1,i2} TCD(i1,i2) D(i1,i2)/abs(D(i1,i2) <\phi_{i1}|\partial H/\partial R|\phi_{i2}>$, here R is the nuclear coordination. Then to conserve the total energy we have used: $V(R)=V(R)+ x F(R)/m(R)$, here m(R) is the nuclear mass, and x is chosen for the smallest value which satisfies the total energy conservation. This x is reported in the screen output of PWmat in the line: "TDDFT,boltzmann,kin:istep,x,dE,a*x**2+b*x+c=0", together with dE(eV), which is the energy drop due to the Boltzmann factor. The a*x**2+b*x+c=0 is the equation used to solve for x.

**temp**: the fixed temperature when $flag\_scale=2$. Note, the temp1,temp2 in the MD\_DETAIL line will not be used for rt-TDDFT simulation (except the temp1 is used to set an initial velocity of the system). We provide this temp here, to distinguish the temp1 for the initial velocity and the temperature you want to use in the Boltzmann factor.

**tau**: the dephasing time (in fs unit) for all the adiabatic state transitions. If tau=-1 (negative), then one needs to provide an file $IN.BOLTZMANN\_TAU$, which specifies the tau for different states as described below.

**IN.BOLTZMANN\_TAU**: the tau input file when tau < 0 in above. This is recommended for the stability of the TDDFT algorithm. It has the following format:

```bash
    tau($m_1$), tau($m_1$+1),tau($m_1$+2),....., tau($m_2$).
```

There should be $m_2-m_1+1$ numbers in a single line, here the $m_1$ and $m_2$ are the first and last adiabatic state index specified in the $TDDFT\_DETAIL$ line. Note, all the numbers are in fs unit. Also note that, for an adiabatic pair i1 and i2 transition, the actual tau for this tau(i1,i2) is sqrt(tau(i1)*tau(i2)). It is recommended that for the few highest adiabatic states in the adiabatic state window [$m_1$,$m_2$], use very small tau. E.g., for the last few tau close to $m_2$, use 0.0001fs. This will make sure these adiabatic states will not have large amplitude $C(i,j)$, which will help the convergence of the TDDFT algorithm.

There are actually more rigorous formulas to calculate tau(i1,i2) from their eigen energies time dependence:

    $$\tau_{i,j}=2 kT  [<|\partial \epsilon_i(t)/\partial t -\partial \epsilon_j(t)/\partial t|^2>_{ave}]^{-1/2} $$

here $\epsilon_i(t)$ is the eigen energy of adiabatic state $\phi_i(t)$ and the $<>_{ave}$ means for a time average. However, in our simulation, we did not calculate this $\tau_{i,j}$, although one can use the above formulat to estimate the decoherent time between a given pair of adiabatic states.

**istep\_start (optional)**: This is an optional value, to indicate from which TDDFT MD step it begins to use the Boltzmann method. If not input, the default value is 1. One can use this parameter to delay the deployment of Boltzmann method, in order to make the algorithm more stable, or for other purposes. One reason is that, for the first 1 or 2 steps, it is possible the TDDFT step itself is not converged. If the SCF TDDFT is not converged, it might appear to be there are large rotations between different eigen states, then that can introduce large Boltzmann correction, and it is wrong, and make can the algorithm unstable. In that case, one can set istep\_start to be 3 or 5, to make the algorithm more stable. Of course, one can also use this to investigate different physics.

**nstep\_CG (optional)**: This is an optional parameter to control the number of steps in an internal conjugate gradient linear solver. Note, to input this parameter, one has to have the istep\_start paremeter before nstep\_CG. The default value of nstep\_CG is 1000. However, for large systems, the solution of a linear equation can significantly slow down the calculation. One can test the use of nstep\_CG=500, or 250 to speed up this step.

### NAMD\_DETAIL

|Tag|NAMD\_DETAIL|
| --- | --- |
|**Format**|NAMD\_DETAIL = $m_1$, $m_2$, nstep\_out|
||NAMD\_DETAIL = $m_1$, $m_2$, nstep\_out, icycle, nstep\_cycle, icrossk, hc|
|**Default**|None|

This line is needed (either the first format, or the second format) when JOB=NAMD.

Note, when JOB=NAMD, besides NAMD\_DETAIL, it also reads in parameters from MD\_DETAIL, and other possible options from TDDFT\_SPACE, TDDFT\_TIME, TDDFT\_STIME, NAMD\_SPECIAL.

In the NAMD calculation, it performs an conventional Born-Oppenheimer MD, but outputs the wave function overlap between consecutive time steps for the adiabatic eigenstates within the window $[m_1,m_2]$. It uses the MD parameters from MD\_DETAIL. The output is written in OUT.NAMD, and will be used in post-process programs like "namd\_dm.x" in \href{http://www.pwmat.com/module-download}{Boltzman-NAMD}, to carry out non-adiabatic MD simulation for carrier dynamics. The carrier wave function $\psi(t)$ will be described by the adiabatic eigen states set within the window $[m_1,m_2]$:

$$
    \psi(t)=\sum_i C_{i}(t)\phi_i(t), i=m1,m2
$$

Note, inside OUT.NAMD (which is an unformatted file), it has the following write-out format:

```fortran
write(10) istep0,islda,nkpt,time,mst_win
do iislda=1,islda
do kpt=1,nkpt
write(10) kpt,iislda,mst_win
write(10) eigen(1:mx,kpt,iislda)
enddo
enddo
do istep=1,nstep
write(10) istep,islda,nkpt,time,mst_win
do iislda=1,islda
do kpt=1,nkpt
write(10) kpt,iislda,mst_win
write(10) hh(1:mst_win,1:mst_win,kpt,iislda)
write(10) eigen(1:mx,kpt,iislda)
enddo
enddo
enddo
```

The first step is different from the rest of the steps. mx is the number of bands calculated,
while $mst\_win=m_2-m_1+1$ is the window of the NAMD output. Eigen is the eigen energy in atomic
unit (Hartree). $hh(m1,m2,kpt,iislda)= <\psi_{m1}^*(istep-1)|\psi_{m2}(istep)>$ for this kpoint and spin iislda.

There are also more advanced options (the second line form). There icrossk=0 or 1: 0, no action for this, 1, output the kpoint cross product in another file: OUT.NAMD\_CROSSK. This is only useful if
multiple kpoints are used. This file has the following format:


```fortran
write(10) istep0,islda,nkpt,time,mst_win
do iislda=1,islda
do kpt=1,nkpt
write(10) kpt,iislda,mst_win
write(10) eigen(1:mx,kpt,iislda)
enddo
enddo
do istep=1,nstep
write(10) istep,islda,nkpt,time,mst_win
do iislda=1,islda
do kpt=1,nkpt
write(10) kpt,iislda,mst_win
write(10) hh2(1:mst_win,1:mst_win,1:nkpt,kpt,iislda)
write(10) eigen(1:mx,kpt,iislda)
enddo
enddo
enddo
```

The only difference is that, in OUT.NAMD, each time step we have: hh(m1,m2), here we have $hh2(m1,m2,kpt2,kpt,iislda)= <u^*(m1,kpt,istep-1)|u'(m2,kpt2,istep>$, so you have cross kpoint dot-product. Note, u is the Bloch part of the wave function, so the cross-k point dot product is not zero. For kpt2.eq.kpt, u' is just the u. But if kpt2.ne.kpt, then u'(istep) has been project out the component of the u(istep-1) in the following way:

    $$u'(m_2,kpt_2,istep)= u(m_2,kpt_2,istep)-\sum_{m_3} h(m_3,m_2,kpt_2)* $$
    $$ [1- exp(-(|h(m_3,m_2,kpt_2)|^4/hc^2)] u(m_3,kpt_2,istep-1)$$

here $h(m_3,m_2,kpt_2)=<u^*(m_3,kpt_2,istep-1)|u(m_2,kpt_2,istep>$. Everything is done within the same iislda.
  
The above output can be used for a special algorithm for postprocess NAMD calculations allowing the hot carrier to jump k-points. One can take $hc$ to be about 0.1, for time step 1fs MD. For better result, one can reduce the time step, while further decrease $hc$.

Another special option is the icycle option. If icycle=0, no action will be taken. If icycle=1, then the program will try to make the MD periodic in time. Note, for this scheme, one cannot restart the MD. So, iMD must equal to 1, not 11, etc. The nstep\_cycle is the buffer region. One possibility, for example, is MDstep=1000 (or 2000), nstep\_cycle=200. The idea is that, the code will make the last step (istep=MDstep) the same (in both atomic position and velocity) as the nstep\_cycle step (note, not the first step). Hence the MD will do a time cycle with a periodicity of MDstep-nstep\_cycle. In the NAMD post-process, this can be used to carry out the NAMD forever, while the nuclear movement has a periodicity of MDstep-nstep\_cycle steps. Note, ath the MDstep-th step, the adiabatic wave function has used the same one as that for the nstep\_cycle-th step. This allows the NAMD calculations to continue forever.

Note, it must be critical to check the potential and kinetic energy in MDSTEPS file, make sure the transition period is smooth. It is critical, for this to work, the whole system from beginning to the end should not drift away. Instead, most atom should only have vibrations (e.g., as in a crystal system).

The nstep\_out is the number of steps interval to output the wave functions (within the window $[m_1,m_2]$) into the ugio.allxxxxxx file. Note, this could be large files, so you probably don't want to output the wave function at every step (e.g., nstep\_out=1). However, sometime the output wave functions can be used to do some postprocessing, which are not done during the NAMD simulation. For example, the wave functions can be used to introduce some additional term $\delta H$, so you have $<\phi_i(t)|\delta H|\phi_j(t)>$ during the postprocess steps. In order to do this, you still don't need to output $\phi_i(t)$ at every time step, because the code does output the overlap $S_{ij}=<\phi_i(t)|\phi_j(t+dt)>$ between consecutive steps, so you can use the $S_{ij}$ to link the adiabatic states at different time. Nevertheless, you might still want to use a relatively small nstep\_out for this regard. As a balance, we found that nstep\_out=50 might be a good choice for many problems. But be prepared for a lot of files!

For more information, please see APPENDIX B.

### NEB\_DETAIL

|Tag|NEB\_DETAIL|
| --- | --- |
|**Format**|NEB\_DETAIL = IMTH, NSTEP, FORCE\_TOL, NIMAGE, AK, TYPE\_SPRING, $E_0$, $E_N$, ITYPE\_AT2, ATOM2.CONFIG|
|**Default**|None|

The NEB\_DETAIL line is needed when JOB=NEB.

For the NEB algorithm, please refer to Ref.\cite{NEB}. In the NEB run, NIMAGE+2 atomic configurations are used, the NIMAGE intermediate configurations connect the initial configuration (in ATOM.CONFIG) with the final configuration (in ATOM2.CONFIG). This is also called the string. During the NEB run, the NIMAGE intermediate configurations will be relaxed together, almost like a NIMAGE*natom atom large system. However, during the atomic relaxation, the atomic force component along the string will be removed (hence not be minimized), so this is call the nudged elastic band method. The goal is to have the force perpendicular to the string to be zero (the string will be moved during the relaxation), while leave alone the atomic force component along the string, meanwhile hopefully keep the distance roughly equal among the NIMAGE+2 image points. Larger the NIMAGE, more difficult is the calculation. For simple problems, typically NIMAGE can be about 5. The output of NEB is written in RELAXSTEPS, and MOVEMENT.

**IMTH**: the algorithm used for atomic relaxation.

|IMTH|Description|
| --- | --- |
|1|conjugate gradient|
|2|BFGS|
|3|steepest decent|
|4|VFF preconditioned conjugate gradient|
|5|Limited-memory BFGS|
|6|FIRE: Fast Inertial Relaxation Engine|

For NEB calculation, IMTH=5 and 6 are recommended for good convergency. Some options can be set for the optimizers in the optional file IN.RELAXOPT(need to set IN.RELAXOPT=T in etot.input). The IN.RELAXOPT is as follows:
```bash
RELAX\_MAXMOVE = 1.0
\verb"    "! max move distance. unit bohr -- for method=1,5,6.
LBFGS\_MEMORY = 30 
\verb"    "! LBFGS storage size -- for method=5. 
FIRE\_DT = 1.0
\verb"    "! initial time step.  unit fs -- for method=6. 
\verb"    "! the max time step for FIRE method is 10*FIRE\_DT. 
RHOWG\_INTER\_TYPE = 1 
\verb"    "! interpolation type for NEB,0--both rho and wave function; 1--rho.
\verb"    "! default = 1, save time by not writing wavefunction to disk 

```

**NSTEP**: the maximum number of line-minimization steps in the relaxation process. This is the NEB steps.

**FORCE\_TOL**: the atomic force tolerance ($eV/{\angstrom}$) to stop the relaxation. This is the maximum atomic force (after the component along the string direction has been projected out) of all the atoms and all the images.

**NIMAGE**: the number of images in the NEB method (these are the images except the initial and final two valleys). So, there are in total NIMAGE+2 configurations in the string of images connection the initial and final configurations. In a NEB calculation, NIMAGE+2 atomic configurations (called images) are used, which connect the configuration from the initial state to the final state. Initially, the NIMAGE intermediate images are generated by linear interpolations of the two end images (the two end images, one initial, one final, are input by the user). The two end images will not be changed, while the NIMAGE intermediate images will be relaxed.

**AK**: the spring constant for the image string ($eV/{\angstrom}^2$). In the NEB, a string connecting the images are used to ensure the coverage between the initial and final configurations. $AK$=0.1 to 1 $eV/{\angstrom}^2$ are reasonable values. Larger $AK$ (especially for TYPE\_SPRING=2), better the convergence, but it can introduce bigger errors (for TYPE\_SPRING=2).

**TYPE\_SPRING**: the type of string used in NEB algorithm.

|TYPE\_SPRING|Description|
| --- | --- |
|1|the original NEB algorithm (where the string force perpendicular to the string tangent is removed)|
|2|a conventional string, the perpendicular string force is not removed|
|3|the regular NEB algorithm(Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points) \cite{RENEB}|

TYPE\_SPRING=2 converges better, but it can introduce an error (larger $AK$, larger the error). But one can first use larger $AK$, then after the initial NEB relaxed, re-runs NEB using smaller $AK$ (or TYPE\_SPRING=1,3). This will help the convergence.

Usually the TYPE\_SPRING=1,2,3 converge good with IMTH=5,6. One can directly try the TYPE\_SPRING=1 or 3 with IMTH=5 or 6, if bad, then follow the above advice(use TYPE\_SPRING=2 first).

Additionally, one can use TYPE\_SPRING=11,22,33 for CI-NEB(11-original CI-NEB, 22-conventional string CI-NEB, 33-regular CI-NEB). The climbing image NEB(CI-NEB) method constitutes a small modification to the NEB method\cite{CINEB}. Information about the shape of the MEP is retained, but a rigorous convergence to a saddle point is also obtained. CI-NEB will choose the image with highest energy as the climbing image, then reformat the force on the climbing image. The force on the climbing image is the full force due to the potential with the component along the elastic string inverted, so the climbing image is not affected by the spring forces, and will climb up along the elastic string to converge rigorously on the highest saddle point. 


Some suggestions

1. How to choose NEB method -- TYPE\_SPRING

As mentioned above, TYPE\_SPRING=2 is easy to converge, but can introduce large error. So if TYPE\_SPRING=1 or 3 converges bad, you can first run the TYPE\_SPRING=2 to get a better guess of the MEP, then rerun NEB with TYPE\_SRPING=1 or 3.

For TYPE\_SPRING=1 (original NEB) and TYPE\_SPRING=3 (regular NEB), both has intrinsic instability. For the original NEB, "In systems where the force along the minimum energy path is large compared to the restoring force perpendicular to the path and when many images of the system are included in the elastic band (string), kinks can develop and prevent the band from converging to the minimum energy path". For the regular NEB, the spring force's formulation omits the the spring force that is perpendicular to the local tangent, then you may encounter that pathway deviate from the MEP and the overall path lengths may grow out of proportion. \cite{AUTONEB}

Maybe your first choice can be TYPE\_SPRING=1, if you find the 'kinks' is the problem, then try TYPE\_SPRING=3.

2. How to choose relaxation method -- IMTH

All IMTH=1,2,3,4,5,6 are methods based on forces, and IMTH=1,2,3,4 each has an exact line-minimization, IMTH=5,6 do not have exact line-minimization. If the forces are not smooth in the calculations, methods with exact line-minimization may not converge well (you will see jumps of the energy as a function of iteration steps). In NEB calculations we recommend the IMTH=5,6. IMTH=5 (Limited-memory BFGS) is fast, and you can set RELAX\_MAXMOVE in file IN.RELAXOPT to make it more stable or more aggressive. IMTH=6 (Fast Inertial Relaxation Engine) is stable, and you can set FIRE\_DT in file IN.RELAXOPT to get more stable or more faster convergence.

We recommend the first choice is IMTH=5.

3. How to choose the pseudopotentials

PWmat provide several types of pseudopotentials, NCPP-SG15, NCPP-PD03, NCPP-PD04, ONCV-PWM-PBE, etc.
We recommend to try to  use the ONCV-PWM-PBE first, it is more smooth, and more easy to converge for relaxation. However, you might want to test this pseudopotential with the more accurate ones, like SG15 and PD03. There are cases the PWM-PBE gives the wrong results. But if it works, you can use it for large system calculations.
If you want to use SG15,PD03,or PD04, we recommend to use bigger Ecut (e.g, 50, 60 Ryd), and set Ecut2=4*Ecut (to avoid egghead problem).

4. How to do CI-NEB

Do not use CI-NEB from the beginning, that will converges bad. PWmat CI-NEB will choose the climbing image automatically, so the climbing image could change during the convergence process if do CI-NEB from the beginning.
We recommend first converge your NEB calculation, then rerun use the CI-NEB. (CI-NEB is specified use TYPE\_SPRING: 11-original CI-NEB, 22-conventional string CI-NEB, 33-regular CI-NEB).

5. The last but not least

Make sure the SCFs are converged.

**E\_0,E\_N**: the precalculated (e.g., using JOB=RELAX) initial ($E_0$) and final ($E_N$) local minima energies (in $eV$) for configurations in ATOM.CONFIG and ATOM2.CONFIG. Actually, these numbers are not used in the algorithm, but will make plotting more straight forward.

**ITYPE\_AT2, ATOM2.CONFIG**: the type of ATOM2.CONFIG file and the atomic position file name: ATOM2.CONFIG.

|ITYPE\_AT2|Description|
| --- | --- |
|1|the type of ATOM2.CONFIG file is the second minimum configuration (the first local minimum configuration is given in IN.ATOM = ATOM.CONFIG). Then, from ATOM.CONFIG to ATOM2.CONFIG, NIMAGE equal distance images will be created by linear interpolations.|
|2|ATOM2.CONFIG contains all the NIMAGE+2 image configurations, including the initial and final images. Thus ITYPE\_AT2=2 is a continued NEB run following the previous NEB runs. You can use the final.config from a previous unconverged NEB run as ATOM2.CONFIG, or copy a group of images from MOVEMENT file. In this case, the ATOM.CONFIG in IN.ATOM = ATOM.CONFIG is not used (but that line still need to be provided). There are cases where the linear interpolation between the first and the last configuration will generate some unphysical NIMAGE IMAGIES (e.g., with atoms too close together, or bond orders wrong). In that case, the user can use ITYPE\_AT2=2 to provide NIMAGE intermediate images manually, to avoid the unphysical interpolation.|

### SCFEP\_DETAIL

|Tag|SCFEP\_DETAIL|
| --- | --- |
|**Format**|SCFEP\_DETAIL = Level1, Level2, $\alpha$, Numkpt, Numspin|
|**Default**|None|

This is a required line for JOB=SCFEP electron-phonon coupling constant calculation.
The calculated electron-phonon coupling constant will be reported in OUT.EP\_COEFF. See JOB=SCFEP (see \ref{jobscfep}) for details. In the JOB=SCFEP calculation, one must use IN.WG=T, an input wave function file will be provided. The electron-phonon coupling constant written in OUT.EP\_COEFF will report: $<\psi(i_1,k,s)|\partial H/\partial R|\psi(i_2,k,s)>$, here $\psi(i,k,s)$ is the input wave function from IN.WG for state index i, kpoint k, and spin s.

**Level1, Level2**: The wave function index $i_1$ and $i_2$

**$\alpha$**: a small number (e.g., 0.1) used to add $\alpha \psi(i1)*\psi(i2)$ onto the charge density to do SCFEP calculation. Here, we assume $\psi(i1)$ and $\psi(i2)$ are real. Suggested $\alpha$ is 0.1. Smaller this number, the numerical derivative will be more accurate, but it also require higher level convergence for the SCF calculations.

**Numkpt**: The kpoint index k for $\psi(i,k,s)$.

**Numspin**: The spin index s (1 or 2) for $\psi(i,k,s)$.

### SCF\_SPECIAL    

|Tag|SCF\_SPECIAL|
| --- | --- |
|**Format**|SCF\_SPECIAL = iflag, Ef0, i1, i2, j1, j2, k1, k2|
|**Default**|SCF\_SPECIAL = 0|

This special option is for nonequilibrium boundary condition calculation for JOB=SCF. Currently, it does not support JOB=NONSCF and JOB=MD (it can run, but the results might not be good). In many problem (for example, device simulation), there could be fixed electrode potentials, there we like the potential to satisfy some specific boundary condition. Note, this might be different from the fixed Fermi Grand canonical calculation, where one electrode potential is fixed. The fixed Fermi is often used together with solvent model with Poisson-Boltzmann method. In the current case, it is often for pure solid calculation, and there is no implicit solvent model. Besides, it is often the case, there are several electrodes. We like the potential on some boundaries (enclosing boundary) to be the given values (e.g., representing the on and off of a CMOS gate, or the source and drain bias potential). Furthermore, this is truely an nonequilibrium simulation, as these boundary values (of the electrode voltage) is often related to the local Fermi energy. So, there is not a single Fermi energy in the simulation. We thus have to use position dependent Fermi energy. Thus, iflag=1 will represent such JOB=SCF calculations. It will give the SCF potential file under such boundary condition. This potential file can be used for quantum transport calculation.

To do such nonequilibrium calculation, we will do two steps.

In the step one, a normal JOB=SCF calculation (iflag=0, or without the SCF\_SPECIAL line) will be carried out. We will take its global Fermi energy as Ef0 (eV) input in the above line. Besides, copy OUT.VR into IN.VR0, which will be used for iflag=1 calculation. Also please copy OUT.RHO into IN.RHO, OUT.WG into IN.WG for
subsequent iflag=1 calculation for fast convergence.

In the second step, set iflag=1. Place Ef0 as mentioned above. Now, we need to use i1,i2,j1,j2,k1,k2 to specify the boundary condition. These are the grid points in the n123 grid of the iflag=0 calculation. More specifically, $i1,i2 \in [1,n_1]$, $j1,j2\in [1,n_2]$, $k1,k2\in [1,n_3]$. These are the planes where the fixted potential will be speficied. Note, the boundary condition is not determined by the edge of the periodic box, instead they are specified by these plane. These plane will define a smaller box inside the periodic box, the corner of this smaller box is at [i1,j1,k1], while the size of this box is: (i2-i1,j2-j1,k2-k1). Note, for this calculation, it is essential to shift the coordinates, so this small box is in the middle of the periodic box. Besides, it is also for some dimention to be periodic. in that case, the corresponding ijk1,ijk2 should both be zero. For example, to have the periodic condition in the second dimension, we should have j1=0,j2=0.

Now, for the dimensions which are not periodic boundary condition (their ijk1,2 are not zero), we need to prepare the corresponding IN.2D\_VR.1, IN.2D\_VR.2, IN.2D\_VR.3 file. For example, if i1,i2 are not zero, then we need a IN.2D\_VR.1 file. This file specifies the dEf(r) on the two planes of i1,i2. More specifically, it is written in the following format (ascii file, so it can be viewed and plotted):

```fortran
    do k=1,n3
    do j=1,n2
    write(IN.2D_VR.1,*) dEf1(j,k), dEf2(j,k)
    enddo
    write(IN.2D_VR.1,*)
    enddo
```
Note, the two columns are for the i1,i2 two planes. There is an empty line, and the for i,j,k (n1,n2,n3), the earlier index are numerated first. This special format is used, so it can be plotted in gnuplot using "splot".

The units for dEf1,dEf2 are in eV. These IN.2D\_VR.1,2,3 files are prepared by utility programs (e.g., gen\_D2V.f), or written by the user. Note, the dEf1,2 are shift of the Fermi energy (also the potential V(r)) in reference to the Ef0 (V0(r)) at those boundary. Thus, the potential V(r) at those boundaries equals to V0(r)+dEf(r). After this, one can run the PWmat again. It will generate the system under the fixed boundary
condition, the OUT.VR can be used to calculate the quantum transport for device simulation.

Here, we explain the underlying algorithm used to carry out the calculation. There are two aspects. One is the Poisson solution to satisfy the fixed boundary condition, another one is the occupation of the wave functions. For the Poisson solution, for a given charge $\rho(r)$, we first use the conventional FFT method to solve a periodic Poisson solution, to get $V_P(r)$. Note, in order to satisfy V(r) equals V0(r)+dEf(r) at the boundary, we can define: $dV_B(r)=V0(r)+dEf(r)-V_P(r)$ on the boundary (of the inner box, defined by i1,i2,j1,j2,k1,k2). Then solve a fixed boundary condition Poisson equation with zero charge, get dV(r) for values inside the box. This is solved by the Fishpack package. Then $V(r)=V_P(r)+dV(r)$. For dV(r) outside the box, we have used simple extension, from its boundary value.

For the occupation, we have used a spatial dependent Fermi energy $Ef(r)=V(r)-V0(r)+Ef0+dE$. The occuption is done with with spatial dependent Fermi energy as: $\rho(r)=\sum_i occ((\epsilon_i-Ef(r))/kT) |\psi_i(r)|^2$, here $\epsilon_i$ is the eigen energy, and occ is the Fermi-Dirac occupation function. Note, we have used a small shift dE to guarantee we get the exact required total charge. Note, this procedure gaurantee the local Fermi energy is at the same position of the local density of state compared with the neutral charge calculation at the step 1. This local Fermi energy is important. For many device systems, the Fermi energy at the electrode is higher than the conduction band in the substrate etc. So, if a global Fermi energy is used, it can be cause large charge slashing from one side to another, which is not physical. In reality, under the open boundary condition (for the current), the system can maintain an steady state nonequilibrium solution, where spatially dependent Fermi energy exists.

### MD\_SPECIAL

|Tag|MD\_SPECIAL|
| --- | --- |
|**Format**|MD\_SPECIAL = iflag, x1, x2, x3, Rcut, dR, dV|
|or|MD\_SPECIAL = iflag, x1, x2, x3, Rcut, dR, dV, frac, P, rate|
|**Default**|MD\_SPECIAL = 0|

**iflag**: the flag for the MD special constraint. See below for details.

|iflag|Description|
| --- | --- |
|0|no special constraint|
|1| this is a special constraint in MD. In this option, at t=0, all atoms within Rcut from the center (x1,x2,x3) will be frozen, and all the other atoms will be subjected to a spherical potential with height dV (eV) centered at (x1,x2,x3) (fractional coordinate) with radius cut-off Rcut (ansgtrom) and a buffer  dR (in a potential file as: dV*exp(-(r-Rcut)/dR)/(exp(-(r-Rcut)/dR)+1)). This is used to keep the atoms out (dV>0) from one domain, or keep the atoms within one domain (dV<0).|
|2|this is the same as iflag=1, except, the atoms within Rcut will not be frozen.|
|22|this is the same as iflag=2, except, when calculating the distance for one atom position (x1a,x2a,x3a) to the spherical center (x1,x2,x3), one do not do periodic wraping for the atoms position. Instead, the atoms position (x1a,x2a,x3a) are defined within ([0,1],[0,1],[0,1]) range. This can be useful, for example, to restraint the water on top of a slab, where one place the (x1,x2,x3) below the slab, so the water will be inside a halfdome on top of the slab, but it will not affect the water on the other side  of the periodic box.|
|3|this option requires the frac,P,rate in the input line. frac is the fraction [0,1] of the halfdome to the whole sphere. E.g., if the sphere is a halfdome above a substate (like in iflag=22,33), then the frac will be less than 1. This is used to estimate the surface of the halfdome. P: the desired pressure in the unit of bar (0.987atmosphere). rate is a MD change rate for each MD step. The idea here is that, during MD, one can change Rcut (radius) of the confinement sphere (usually, dV is negative for this), so the pressure on the confinement wall equals to P. This will help to decide what Rcut one should use. Note, P should not be too small, otherwise it might be difficult to converge. We recommend(e.g.): P = 50 bar (atmosphere), rate=0.02. Note, one should have a reasonably large dR, so the force will not be too large. For example, dV=-2 (eV), dR=1(A).  frac is determined from the geometry. If it is a full sphere constraint, then frac=1. For iflag=33, e.g., an halfdome constraint, then it is the solid angle ratio between the halfdome solid angle and the 4$\pi$ full sphere solid angle (this is used to calculate the halfdome surface area). The actually pressure (P\_sph), the input desired pressure P\_MD\_sp (=P), and the dynamically adjusted Rcut (Rcut\_MD\_sp) will be reported in the header of MOVEMENT file for each step.|
|33|same as iflag=3, except, when calculating the distance for one atom position (x1a,x2a,x3a) to the spherical center (x1,x2,x3), one do not do periodic wraping for the atoms position. Instead, the atoms position (x1a,x2a,x3a) are defined within ([0,1],[0,1],[0,1]) range. This can be useful, for example, to restraint the water on top of a slab, where one place the (x1,x2,x3) below the slab, so the water will be inside a halfdome on top of the slab, but it will not affect the water on the other side  of the periodic box.|

In iflag=1,2,3,22,33, one can use a "Weight\_Atom" section in xatom.config file to specify which atom can feel this potential (or how much weight to feel this potential). More specifically, in the xatom.config file, one can add a section like:

```bash
Weight\_Atom
1
iatom(1), weight(1)
iatom(2), weight(2)
....
iatom(natom),weight(natom)
```

The weight(i), for i=1,natom can be used for many special purposes in the code where a coefficient for an atom is needed. iatom(i) is the atom z-number. So, as a result, the potential felt by each atom is dV*weight(i). So, one can turn-off some of the atoms to feel this particular potential.

### MD\_SPECIAL2

|Tag|MD\_SPECIAL2|
| --- | --- |
|**Format**|MD\_SPECIAL2 = iflag, Rcut, dR, dV|
|**Default**|MD\_SPECIAL2 = 0|

This is typically used together with MD\_SPECIAL=1, with its center defined using x1,x2,x3 defined in the MD\_SPECIAL line. It defined another potential using parameters Rcut, dR and dV. However, it will be applied to the atoms with weight(i,2), which is defined in the xatom.config file section Weight\_Atom as:

```bash
Weight\_Atom 
2
iatom(1), weight(1), weight2(1)
iatom(2), weight(2), weight2(2)
....
iatom(natom),weight(natom), weight2(natom)
 
```

### MD\_SPECIAL3

|Tag|MD\_SPECIAL3|
| --- | --- |
|**Format**|MD\_SPECIAL3 = iflag, m1, m2|
|**Default**|MD\_SPECIAL3 = 0|

This is a very special case for electron phonon coupling calculation. The idea is to input m1 wave function $\phi_1(i)$ for i=1,m1, and m2 wave function $\phi_2(i)$ for i=1,m2 at the beginning of the MD simulation from files IN.WG1\_MDSP3, IN.WG2\_MDSP3. Then during the MD simulation, at every dt step, it will output: $hh(i,j)=<\phi_1(i)|H(t)|\phi_2(j)>$. Note, the H(t) is Hamiltonian at time t, and $\phi_1(i)$ and $\phi_2(j)$ are not changed during the MD.

The $m1\times m2$ double complex matrix $hh(i,j)$ is written inside the binary file: MDSP3.hh.out 
in an concatenation style (position="append"), i.e, continuesly written at its end. It has
the following format:

```bash
t, m1, m2, nkpt, islda
hh(m1,m2,nkpt,islda)
t, m1, m2, nkpt, islda
hh(m1,m2,nkpt,islda)
....
t, m1, m2, nkpt, islda
hh(m1,m2,nkpt,islda)
```

One can used the small utility file plot\_MDSP3.f to plot it out. The IN.WG1\_MDSP3, IN.WG2\_MDSP3 have the same format as IN.WG (if there are spin=2, then we also need IN.WG1\_MDSP3\_2, IN.WG2\_MDSP3\_2). But their number of wave functions m1, m2 can be much smaller than mx in IN.WG. They can be selected from OUT.WG from some previous runs (e.g, defect wave functions) using utility file: Select\_WG.f90. 

### MD\_SPECIAL4

|Tag|MD\_SPECIAL4|
| --- | --- |
|**Format**|MD\_SPECIAL4 = iflag, Rc|
|**Default**|MD\_SPECIAL4 = 0|

This is also usually used together with MD\_SPECIAL=1, but with solvent model. It excludes the solvent effect from the center of the x1,x2,x3 defined in the MD\_DPECIAL=1 line, with a radius of Rc (A unit). This is done by adding some charges at the center when defining the dielectric constant. This is used to prevent the solvent effects to intrude into the center region.

### NAMD\_SPERICAL

|Tag|NAMD\_SPERICAL|
| --- | --- |
|**Format**|NAMD\_SPERICAL = iflag\_NAMD\_sp, param1, param2|
|**Default**|NAMD\_SPERICAL = 0|

This is an optional input section for special input parameters (param1,param2) for JOB=NAMD calculations. Default iflag\_NAMD\_sp=0, the parameters are not used.

When iflag\_NAMD\_sp=1, it will output the dipole moment matrix at every MD step in a file called OUT.NAMD\_SP. It is a binary file, with the following format.

```fortran
Do istep=1,nstep

write(OUT.NAMD\_SP) istep,islda,nkpt,param1,param2,mst\_win,time

Do iislda=1,islda

DO kpt=1,nkpt

write(OUT.NAMD\_SP) PXYZ

ENDDO

ENDDO

ENDDO

```

Here PXYZ is a complex*16 matrix PXYZ(mst\_win,param2-param1+1,3). mst\_win=m\_2-m\_1+1 from NAMD\_DETAIL, are the number of states output for NAMD calculation. Thus: $PXYZ(i,j,k)=<\phi_i|P_k|\phi_j>$, here k=1,2,3, and $P_k$ is the momentum operator $i\nabla_k$.

### CHARGE\_DECOMP

|Tag|CHARGE\_DECOMP|
| --- | --- |
|**Format**|CHARGE\_DECOMP = T / F|
|**Default**|CHARGE\_DECOMP = F|

This will create an atomic charge on each atom using the Hirshfield algorithm, and reported the result in output file OUT.QDIV for JOB = SCF calculation. For spin = 2, it will report total charge, and magnetic\_moment = charge\_up - charge\_down. If the job is to do molecular dynamics (JOB = MD), then the atomic charge will also be reported in MOVEMENT. This option is useful for charge analysis. For the Hirshfield algorithm, the charge on one atom i is defined as: $Q_i=\int \rho(r) {\rho_{atomtype(i)}(r-R_i) \over \sum_j \rho_{atomtype(j)}(r-R_j)} d^3r $, here $\rho_{atomtype(i)}(r)$ is the neutral charge density of the atom type atomtype(i) described in the atom pseudopotential file xxx.upf. Note, this option only works for norm conserving pseudopotential.

### ENERGY\_DECOMP

|Tag|ENERGY\_DECOMP|
| --- | --- |
|**Format**|ENERGY\_DECOMP = T / F|
|or|ENERGY\_DECOMP = T / F, type|
|**Default**|ENERGY\_DECOMP = F 1|

This will decompose the total DFT (LDA, PBE only, but it also works for IN.SOLVENT=T, as well as Poisson-Boltzmann equation) energies into the energies belong to each atom (atomic energies). The sum of the atomic energies will be equal to the total DFT energy (but differ by a constant, this constant is independent of the position of the atom, as well as the lattice length. Basically, each atom will miss an atom type specific energy constant, an onsite energy term). It rewrites the total energy as an spatial integral of the positive energy density term, and use the Hirshfield algorithm to decompose such energy density into each atom, much like the above decomposition for the charge density. See Ref.[J. Kang, L.W. Wang, Phys. Rev. B 96, 020302(R)(2017)]  for details. For JOB = SCF, the decomposed energy will be reported in OUT.ENDIV. For JOB = MD, the decomposed energy will also be reported in MOVEMENT. These decomposed energies can be used to do force field fitting.

**type**: This is an optional input. The default is type=1. There are four options for type,
|type|Description|
| --- | --- |
|1,11|the electrostatic energy density is expressed as: $1/8\pi \|E(r)\|^2$ (here E is the electric field). This is positive everywhere, but it can have large amplitude even in vacuum region.|
|2,22|the electrostatic energy density is expressed as: $1/2 \rho(r) V(r)$, here $\rho(r)$ is the charge density including both electron and nuclear charge, V(r) is the total electrostatic potential. So, this only has values where $\rho(r)$ is not zero. Note, for all types, for IN.SOLVENT=T, the solvent polarization induced electrostatic energy density is always represented as $\rho_{solute}(r)V_{polarization}(r)$.|
|1,2|a straight forward Hirshfeld spatial partitioning is used to partition the energy density, to yield the energy for each atom. However, for the Hirsheld partitioning, the atomic charge (amplitude and shape) can be altered by the parameter in ENERGY\_DECOMP\_SPECIAL and ENERGY\_DECOMP\_SPECIAL2.|
|11,22|an atomic weight watom(iatom) is used, and dynamically adjusted, so when doing charge decomposition, it will yield the atomic charge equal to the neutral atom charge z(atom\_type). This watom(iatom) is then used to do the energy decomposition. The hope is that, by getting the fixed charge density, the energy part (especially when using type=22) will have minimum variations. All these are designed to get the minimum fluctuation of the atomic energy. Note, there are additional costs by doing type=11,22, in order to find watom(iatom) for each atomic configuration. The watom(iatom) are listed in the last column in OUT.ENDIV, or the corresponding section in MOVEMENT.|

>
>**TIP**: If ENERGY\_DECOMP\_COULOMB=T, then type=1/2 are the same, and type=11/22 are the same. The detailed option will be determined by imth in the ENERGY\_DECOMP\_COULOMB line. Please check that section.
{: .block-tip}

>
>**WARNING**: ernergy decomposition does not support NUM\_BLOCKED\_PSI = T .
{: .block-warning}

### ENERGY\_DECOMP\_SPECIAL

|Tag|ENERGY\_DECOMP\_SPECIAL|
| --- | --- |
|**Format**|ENERGY\_DECOMP\_SPECIAL = w(1), w(2), ...., w(ntype)|
|**Default**|ENERGY\_DECOMP\_SPECIAL = 1, 1, ..., 1|

This is an additional option for ENERGY\_DECOMP=T as well as CHARGE\_DECOMP=T. This option modifies the weight for each atom type, not just using $\rho_{atomtype(i)}(r)$, but using $w(atomtype(i))*\rho_{atomtype(i)}(r)$ as the weight in the Hirshfeld method to calculate the charge and energy. The default values for all w(i) are 1.

### ENERGY\_DECOMP\_SPECIAL1

|Tag|ENERGY\_DECOMP\_SPECIAL1|
| --- | --- |
|**Format**|ENERGY\_DECOMP\_SPECIAL1 = eta(1), eta(2), ...., eta(ntype)|
|**Default**|ENERGY\_DECOMP\_SPECIAL1 = 1, 1, ..., 1|

This is a very special input, please don't use it if you don't know the technical detail.  In order to do the energy decomposition, for each input atom.UPF pseudopotential file, it will generate the corresponding: atom.UPF.ionrhoR, atom.UPF.rhoq, atom.UPF.rhoatom, files. They are the real space $v\_{loc}(r)$ file before and after fitting; q-space $v_{loc}(q)$ file before and after fitting, and $\rho(r)$ used for spatial partitioning function to generate the atomic quantities. It is worth to plot $v_{loc}(r)$, especially $v_{loc}(q)$. In $v_{loc}(q)$, for $0 < q < q_{c}/2$, it is the original $v_{loc}(q)$, while for $q_{c}/2 < q < q_{c}$ is the fitted one. Make sure the fitted one is close to the original one, there is no big variation. If there are large change, one can modify eta(itype). The default eta value is 1. Larger eta can make the $v_{loc}(q)$ smoother, but could be less accurate in some sense. One can set eta to 1.5 for example if $v_{loc}(q)$ for $q_{c}/2 < q < q_{c}$ is large.

This is an additional option for ENERGY\_DECOMP=T as well as CHARGE\_DECOMP=T. This option modifies the weight for each atom type, not just using $\rho_{atomtype(i)}(r)$, but using $w(atomtype(i))*\rho_{atomtype(i)}(r)$ as the weight in the Hirshfeld method to calculate the charge and energy. The default values for all w(i) are 1.

### ENERGY\_DECOMP\_SPECIAL2

|Tag|ENERGY\_DECOMP\_SPECIAL2|
| --- | --- |
|**Format**|ENERGY\_DECOMP\_SPECIAL2 = exp\_decomp a(1) a(2) ... a(ntype) b(1) b(2) ... b(ntype)|
|**Default**|ENERGY\_DECOMP\_SPECIAL2 = 1|

This is used to control the atomic charge density $\rho_{atom}(|r-R|)$ to be used in the Hirshfeld formula.

$\rho_{atom}(|r-R|)=\rho_{atom}(|r-R|) + a(atomtype(R)*exp(-(|r-R|/b(atomtype(R)))^2))$

The together with the weight input in ENERGY\_DECOMP\_SPECIAL, the Hirshfeld partition function for atom R at position r, is given as:

$\rho_{atom}^{exp\_decomp}(|r-R|) w(atomtype(R))/\sum_{R'} \rho_{atom}^{exp\_decomp}(|r-R'|)w(atomtype(R'))$

(Note, there could be an additional weight watom(iatom) for type=11,22 partitioning method). So, higher value of exp\_decomp (e.g., 2) can make the partitioning more local, and the interface more abrupt, but it can also be less smooth. 

Note, the parameters in ENERGY\_DECOMP\_SPECIAL1, and ENERGY\_DECOMP\_SPECIAL2 can also be used to control the CHARGE\_DECOMP.

### ENERGY\_DECOMP\_COULOMB

|Tag|ENERGY\_DECOMP\_COULOMB|
| --- | --- |
|**Format**|ENERGY\_DECOMP\_COULOMB = T / F, iconstr,imth,fac1,fac2,numG,q2W,iout|
|**Default**|ENERGY\_DECOMP\_COULOMB = F|

This is an option for Coulomb potential charge fitting. It can either do it on the flight (imth=1,2), or use an already fitted spherical atomic charge model (funcq\_atom.fit), and recalculate the electrostatic Coulomb interaction energy. This will affect the electrostatic interaction involving both the electron charge and nuclear charge. The idea is to fit the electron charge density, and calculate the $rho_{fit}*rho_{fit}$ interactions analytically (using atom-center pair interaction), and we will have a residuation charge: $rho_{res}=rho-rho_{fit}$, and hopefully this residue is small, and it will not have long range interaction. In terms of fit, we can do on the flight, as in imth=1,2 (this is done for using 1 or 2 Gaussians for each atom, and either we have fitting, or we have no fitting, just use the previously fitted results), or we can have an more extensive fitting, using the vion\_coulomb\_fit.f90 utility file, to generate the funcq\_atom.fit files to be read by the program.

>
>**TIP**: Note, ENERGY\_DECOMP\_COULOMB = T does not work with Solvent model.
{: .block-tip}

**imth=1,2**: fitting the charge with Gaussian on the flight. The idea is to use one or two (numG) Gaussians to represent a charge density at a given atom. **iconstr=1**: (this is only used under imth=1,2), one must provide a IN.CONSTRAINT\_COULOMB to give information for constrains during the density fitting. This can be used to specified, the sum of the charge of a few atoms (or one atom) must be a given number. The fitted charge will output in: OUT.natom\_coulomb\_fit.

When ENERGY\_DECOMP\_COULOMB = T, one "ENERGY\_NATOM\_COULOMB" section must be provided in atom.config file. This section has the following format (note for imth=3, this section is not really used, but nevertheless, please provide an place holder faked section ):

```bash
ENERGY_NATOM_COULOMB
 4                       # natom_fit
 151  15  0.20 0.5       # atom order in xatom; zatom, a1, a2 (A)
 152   9  0.15 0.4       # atom order in xatom; zatom, a1, a2 (A)
 153   9  0.15 0.4
 154   9  0.15 0.4
 ------------------------
 150                                 # natom_fix
 1   8  0.3 0.5 -10.674, 11.914      # atom_order,zatom,a1,a2(A),Q1,Q2
 2   8  0.3 0.5 -10.674, 11.914      # atom_order,zatom,a1,a2(A),Q1,Q2
 3   8  0.3 0.5 -10.674, 11.914      # atom_order,zatom,a1,a2(A),Q1,Q2
 4   8  0.3 0.5 -10.674, 11.914      # atom_order,zatom,a1,a2(A),Q1,Q2
 ..... (150 lines)
```

Note, the a1,a2 are the size of the two Gaussian (even if numG=1, only one Gaussian is used, please provide two columns, for the place holder, the format is fixed). The atom\_order is the index of the atom in the original xatom.config (Note, if the atoms are not in consequetive order for different atom types, this first index will be changed, to the index in output xatom.config file. The natom\_fit is the atoms to be fitted (to get their charge parameters Q1,Q2), while the natom\_fix is the atoms which already fitted before, thus already know the Q1,Q2 (the charge on these two Guassian functions).

natom\_fit and/or natom\_fix can be zero (for imth=3, one can set both to be zero). Note, if desired, it is okay for only fiting the charge for a few atoms, instead of all the atoms. Or one can first fit the Q1,Q2 from some other systems, then use them as natom\_fix, and fit some additional atoms in this system.

**numG=1, or 2**: In above, in the section of ENERGY\_NATION\_COULOMB, there are always two Gaussians. Actually, one can adjust the number of Gaussian, numG=1 means only use 1 Gaussian (but there should still have two columns in the ENERGY\_NATION\_COULOMB section, only the first column is used). When numG=2, two Gaussians are used.

**iconstr=1**: the IN.CONSTRAINT\_COULOMB section has the following format:

```bash
2                       # number of constrains (followed by num lines)
3  1.0 1,2,3            # Num_atom; Qtot, ind1,ind2,ind3,...
1 -1.0 4                # Num_atom; Qtot, ind1,ind2,ind3...
```

Note, each line is one constrant.  Qtot is the total charge of these few atoms within this constraint. Ind1,ind2,ind3 are the atom number index of the atoms within the natom\_fit sequence in the ENERGY\_DECOMP\_COULOMB section of the xatom.config file. Note, they are not the index in the original xatom.config atom sequence. They are index within the list of natom\_fit (e.g., must be less or equal to natom\_fit. For example, the 1, 2, 3 atoms in the above example correspond to atoms 151, 152, 153.

The fitted natom\_fit Q1,Q2 results are shown in OUT.natom\_coulomb\_fit  (together with the natom\_fix Q1,Q2 results input within the ENERGY\_NATOM\_COULOMB section of xatom.config).

**imth=1**: In this method, $\rho_{res}(r)=\rho(r)-\rho_{fit}(r)$, note, $\rho(r)=\rho_{el}-\rho_{nuclear}$. So, the fitting is done for the total charge density (including nuclear), not just for the electron charge density. Then: Electrostatic potential density equals: $E_{elst}(r)={1\over 8\pi} |\nabla V_{res}(r)|^2+\rho_{fit}V_{res}(r)$ (excluding $\rho_{fit}*\rho_{fit}$ interaction). In the Hirshfeld parititiong, only the first term is partitioned, the second term is directly integrated for each atom, thus we have a $E(Q*V_{res})$ term (Q is the fitted charge), which is listed as one column in OUT.ENDIV.

To get the atom decomposed energy  without new fitting(after excluding the fitted charge Coulomb interactions), one can set natom\_fit=0 in the ENERGY\_NATOM\_COULOMB section of xatom.config file.

**imth=2**: In this option, the energy density (excluding the $\rho_{fit}*\rho_{fit}$ interaction) is expressed as: ${1\over 8\pi} [|\nabla V(r)|^2 fac1- |\nabla V_{fit}(r)|^2 fac2]$. Here V(r) is the Coulomb potential of $\rho(r)=\rho_{el}-\rho_{nuclear}$, and $V_{fit}$ is the Coulomb potential of $\rho_{fit}$. So, normally, **fac1 and fac2** should be one. They are provided here, just for the purpose of analysis. Note, in this option, there is no $E(Q*V_{res})$ term, and this column in OUT.ENDIV will be zero.

**imth=3**: In this option, there is no on-the-flight fitting. Instead, a prior fitted atomic charge will be used. In this option, all the atoms need to be fitted, not just for a selected subset. The fitted spherical charge density for each atom type in q-space is input from a file called: funcq\_atom.fit. This funcq\_atom.fit is obtained from the utility code vion\_coulomb\_fit.f90, based on the output charge density file OUT.rho\_EpN ($\rho_{el}-\rho_{nuclear}$)  from a previous ENERGY\_DECOMP\_COULOMB run. The electrostatic energy denisty is expressed as ${1\over 8\pi} |\nabla V_{res}(r)|^2+\rho_{fit}V_{res}(r)$. But unlike in imth=1, the second term is not integrated for each atom, instead, it is included in the Hirshfeld partitioning. This imth=3 can significantly reduce the total electrostatic interaction energy (e.g, by a factor of 1000 in the case of melted NaCl).

**q2w**: a parameter (unit 1/A) used for a factor $exp(-(q*q2w)^2/4)$, this factor is used in the on-the-flight fitting of the Gaussian charge density to reduce the electrostatic energy: $\sum_q |\rho_{res}(q)|^2*exp(-(q*q2w)^2/4)*4\pi/q^2$. The idea is that, we can use this parameter to emphasize only the small reciprocal vector q components, hence only the long range part of the Coulomb interaction. Note, this factor is also used when generating output "OUT.Coulomb\_EpN" and "OUT.Coulomb\_residual". It can be used to filter out the high energy components.

**iout= 0 or 1**: 1 will mean there will be output (every MD step, rewrite) charge density: OUT.rho\_EpN  ($\rho_{el}-\rho_{nuclear}$) on a double grid; OUT.rho\_EpN.fit ($\rho_{fit}$); OUT.Coulomb\_EpN (the Coulomb potential of $\rho_{el}-\rho_{nuclear}$, subject to the $exp(-(q*q2w)^2/4)$ factor; OUT.Coulomb\_residual (the COulomb potential of $rho_{res}=(\rho_{el}-\rho_{nuclear})-\rho_{fit}$ subject to the $exp(-(q*q2w)^2/4)$ factor. Note, these quantitites are double grid ($2n_1\times 2n_2\times 2n_3$) quantities for Hirshfeld partitioning usage. iout=1 can be expensive.

It usually only used for JOB=SCF calculation, and it is can be used to generate OUT.rho\_EpN to fit the funcq\_atom.fit using vion\_coulomb\_fit.f90.


### ATOMIC\_ORBITAL\_IATOM\_OUT

|Tag|ATOMIC\_ORBITAL\_IATOM\_OUT|
| --- | --- |
|**Format**|ATOMIC\_ORBITAL\_IATOM\_OUT = atomic\_orb\_iatom\_i, atomic\_orb\_iatom\_e|
|**Default**|NO DEFAULT|

ATOMIC\_ORBITAL\_IATOM\_OUT contains two parameters atomic\_orb\_iatom\_i and atomic\_orb\_iatom\_e which decide the atom index range(the index is from the IN.ATOM configuration file) for JOB=ATOMIC\_ORB.

## Corrections \& constraints tags

### VDW

|Tag|VDW|
| --- | --- |
|**Format**|VDW = NONE / DFT-D2 / DFT-D3 /PAIR|
|**Default**|VDW = NONE|

This parameter is used to specify the type of Van Der Waals correction.

If use DFT-D2, some variables are optional to be set: LONDON\_S6, LONDON\_C6, LONDON\_RCUT. We use the Grimme’s empirical vdw functional term. It is okay without setting the LONDON parameters.

LONDON\_S6: Global scaling parameter for DFT-D. Default is 0.75.

LONDON\_C6: It is an array its dimension is the number of atomic type. Its format is like this: LONDON\_C6(1) = ..., LONDON\_C6(2) = ... (1),(2) are the atom types, in accordance with IN.PSP1, IN.PSP2. These are the C6 parameters in the Lennard-Jones potential $1/r^6$ term parameter. Only the attractive $1/r^6$ term is included in VDW. The repulsion part is already in the DFT energy. The parameter LONDON\_S6 determines trunction of this term at small r.

The default value is from the Grimme-D2 values. You can refer to the article: S. Grimme, J. Comp. Chem. 27, 1787(2006).

LONDON\_RCUT: The cutoff radius (a.u.) for dispersion interactions calculations. The default is 200.0. For $\left|R1-R2\right|$ larger than this cut-off, the vdW interaction will not be calculated. Note, this default value might be too large.

If use DFT-D3(zero-damping method), some variables are optional to be set: DFTD3\_S6, DFTD3\_RS6, DFTD3\_S18, DFTD3\_RS18, DFTD3\_ALPHA6, DFTD3\_VERSION, DFTD3\_3BODY. For more information about the DFT-D3, one can refer to this paper \cite{dftd3-zero-damping}.

In the D3 correction method, the following vdW-energy expression is used:

$E_{disp} = -\frac{1}{2} \sum\limits_{i=1}^{Natom} \sum\limits_{j=1}^{Natom} \sum\limits_{L}\bigg(f_{d,6}(r_{ij,L})\frac{C_{6ij}}{r_{ij,L}^{6}} + f_{d,8}(R_{ij,L})\frac{C_{8ij}}{r_{ij,L}^{8}} \bigg)$

The dispersion coefficients $C_{6ij}$ are geometry-dependent as they are adjusted on the basis of local geometry (coordination number) around atoms i and j. In the zero-damping method, damping of the following form is used:

$f_{d,n}(r_{ij}) = \frac{s_{n}}{1+6(r_{ij}/(s_{R,n}R_{oij}))^{-\alpha_{n}}}$

where $R_{oij} = \sqrt{\frac{C_{8ij}}{C_{6ij}}}$, the parameters $s_{R,8}$ are normally 1. Respectively, the $s_{6}, s_{8}, s_{R,6}, \alpha_{6}$ are adjustable parameters whose values depend on the choice of exchange-correlation functional. Note the default parameter is tested for PBE.

DFTD3\_VERSION The parameter for zero-damping DFT-D3 method is 3.

DFTD3\_3BODY The parameter controlling whether considering the three body term in DFT-D3 correction method. The default is T (considering). If not considering the term, setting DFTD3\_3BODY = F.

DFTD3\_CUTOFF The cutoff radius (a.u.) for pair interactions in real space. Default is 94.868 bohr.

DFTD3\_CUTOFF\_CN The cutoff radius (a.u.) for coordination number in real space. Default is 40.0 bohr.

VDW=PAIR In that case, the file IN.VDWPAIR need to provided. This is a general user provided pair potential in a numerical form, which has the following format:

```bash
1001, 3, 8.0       ! nr, n_pair, r_cut(angstrom)
31 31 33           ! iatom_1 (npair)
33 31 33           ! iatom_2 (npair)
0.00,  p1,  p2,  p3  ! r, pot(pair1), pot(pair2), pot(pair3) (eV)
...
r,   p1, p2, p3       ! The nr-th line for p1,p2,p3
```

### COULOMB

|Tag|COULOMB|
| --- | --- |
|**Format**|COULOMB = 0 / 1, X1, X2, X3 / 11, X1 / 12, X2 / 13, X3|
|**Default**|COULOMB = 0|

Control the Poisson equation solution (for the Coulomb interaction).

We provide special ways to calculate the Coulomb potential (also called Hartree potential) of the charge density $\rho(r)$ during the SCF calculation. This could be particularly helpful for isolated system calculation, or for slab calculation. This is to avoid electrostatic image interaction for isolated systems, or the dipole moment effect for neutral slab calculations.

|COULOMB|Description|
| --- | --- |
|0|the periodic boundary condition, the default.|
|1, X1, X2, X3|the isolated cluster boundary condition. It can avoid the image interaction in this calculation. The X1, X2, X3 (value: 0~1) are the fractional coordination values in the unit cell edge vectors 1, 2, 3, used to cut a box for this special Coulomb solution. In the other word, the center of the box is at: (X1+0.5, X2+0.5, X3+0.5).|
|11, X1|A slab calculation along the first direction, with the cut at X1. This can avoid the dipole moment interaction between slabs. Note, this method only works for neutral system. For charged slab system, the energy is infinite for an isolated slab. In that case, please just use COULOMB=0.|
|12, X2|A slab calculation along the second direction with the cut at X2. Only for neutral system.|
|13, X3|A slab calculation along the third direction with the cut at X3. Only for neutral system.|

### LDAU\_PSP

|Tag|LDAU\_PSP|
| --- | --- |
|**Format**|LDAU\_PSP1 = LDAU\_L(1), Hubbard\_U(1) (optional Hubbard\_U2(1))|
||LDAU\_PSP2 = LDAU\_L(2), Hubbard\_U(2) (optional Hubbard\_U2(2))|
||...|
|**Default**|LDAU\_PSP1 = -1|
||LDAU\_PSP2 = -1|
||...|

If this parameter is set, LDA+U method will be used. When using LDA+U method, one must specify, for each element(i), the atomic orbit to add U, and the value of U.

Note the (i) should correspond to the IN.PSP(i) for the pseudopotential input.

**LDAU\_L(i) = -1/0/1/2/3**, 0/1/2/3 means adding a U term to the s/p/d/f orbital. -1 means not to use LDA+U.

**HUBBARD\_U(i)**: the U parameter ($eV$) for species element types i, the default value is 0.0. If HUBBARD\_U2(i) is also provided, then the first number is for spin up component, the second number is for spin down component. They can be different. If the second number is not provided, then the spin up and down U parameters will be the same.

Note, there are cases where the same atom type, say Co, needs to use different U parameter depending on its local environment and valence state (e.g., Co$^{3+}$ and Co$^{2+}$). In such case, one can change the atom number of Co in the atom.config file, say, one is 27, another is 127. Also, one add another Co pseudopotential Co2.xxx.upf, and in Co2.xxx.upf, change its atomic number from 27 to 127. Now, one can add different U for this new 127 Co type.

Note, the LDA+U calculation is done by the following Hamiltonian:

$$ H = H_{LDA} + \sum_{I,\sigma} {U_{I,\sigma}\over 2} Tr[n^{I,\sigma}(I-n^{I,\sigma})] $$

here I denote the atomic site, and $\sigma$ is for the spin, and $n^{I,\sigma}(m1,m2)$ is an occupation matrix for atomic orbital $\phi_m$ as:

$$n^{I,\sigma}(m1,m2) = \sum_j <\psi_{j,\sigma} |\phi^I_{m1}> < \phi^I_{m2}|\psi_{j,\sigma}> occ(j,\sigma) $$

Here $\psi_{j,\sigma}$ is the Kohn-Sham orbital for spin $\sigma$ and $occ(j,\sigma)$ is its occupation. In above $U_{I,\sigma}$ is provided by the LDAU\_PSP1 lines.

However, we can also add another term, which is:

$$ H = H_{LDA} + \sum_{I,\sigma} {U_{I,\sigma}\over 2} Tr[n^{I,\sigma}(I-n^{I,\sigma})] + \sum_{I,\sigma} \lambda_{I,\sigma} Tr[n^{L,\sigma}] $$

The parameters $\lambda_{I,\sigma}$ are provided in a special section LDAU\_lambda in the atom.config file:

```bash
LDAU_lambda
    27  0.5  0.5   ! iatom, lambda_up, lambda_dn
    27  0.0  0.0
    ...........
    27  0.0  0.0   ! must have natom lines
```

Note, in this LDAU\_lambda section in atom.config, even if not all atoms have LDA+U, you must provide natom (all atom) lines for every atom. For those atoms which do not has LDA+U (determined by the LDAU\_PSPx lines), their corresponding lambda will not be used. So, if you want to use lambda, you must set the LDAU\_PSPx line for that atom type. You can set a very small U parameter for that atom type to reduce the first LDAU term in above H formula.

In the LDAU calculation, we will have output: OUT.LDAU\_NS, which writes out the $n^{I,\sigma}$ matrix for each atoms which has the LDAU\_PSPx line. One can use this $n^{I,\sigma}$ and $\lambda_{I,\sigma}$ to calculate U with the linear-response method (please check the corresponding PWmat module).


### LDAU\_RCUT\_PSP

|Tag|LDAU\_RCUT\_PSP|
| --- | --- |
|**Format**|LDAU\_RCUT\_PSP1= rcutu1|
||LDAU\_RCUT\_PSP2= rcutu2|
||...|

The Rcut of each element type for the LDA+U calculation for the nonlocal projector. The unit is in Bohr, not Amgstron. This is like the IN.PSP\_RCUTi, but instead of for the full atomic wave function, it is for the nonlocal projector. Note, usually, these rcutui should be bigger than rcuti. If no explicit input for LDAU\_RCUT\_PSPi is provided, the default value of 6 Bohr is used (which is rather large). Note, for LDA+U calculation, it might be necessary to use a large rcutui, e.g., 6, to get a fully converged and smooth force (e.g., for better RELAX convergence).  However, if Ecut2 is very large, then a slightly smaller rcuti can be used.

### STRESS\_CORR

|Tag|STRESS\_CORR|
| --- | --- |
|**Format**|STRESS\_CORR = $num\_pw_1,energy_1,num\_pw_2,energy_2$|
|**Default**|NO DEFAULT|

As we know, the number of plane wavefunction has respond to the size of the lattice, i.e. different lattice will give different number of plane wavefunctions in the same cutoff. When running cell relaxation, the lattice will change. Correspondingly, the number of plane wavefunctions should also change. However, the change of the plane wavefunctions will disturb the procedure of the relaxation. As a result, in DFT calculations, we keep the number of plane wavefunctions all the same. In order to overcome the problem in cell relaxation, we implement an stress correction. One can refer to the paper \cite{stress1} for more details.

The steps to carry out the stress correction: Before running the cell relaxation, you should do two SCF calculations with different cutoff. Each calculation will give the $num\_pw$ (See ``Weighted average num\_of\_PW for all kpoint'' in REPORT), the $energy$ (See ``E\_tot'' in REPORT). Then continue doing the cell relaxation with these parameters setting STRESS\_CORR.

>
>**WARNING**: The two cutoffs should be close enough(within difference of 1~2 ryd ), otherwise the stress correction will not work.
{: .block-warning}

### FIX\_FERMI

|Tag|FIX\_FERMI|
| --- | --- |
|**Format**|FIX\_FERMI = T/F, E\_Fermi, mix\_Q, drho\_pulay|
|**Default**|FIX\_FERMI = F|

This control indicates whether to use fix Fermi energy (fix electrode potential) calculation. The default is F (use fixed total charge). Note, for T, this is usually used together with IN.SOLVENT=T, and with POISSON\_BOLTZMANN = T (for Poisson-Boltzmann screening) in the IN.SOLVENT file for electrochemistry grand cannonical calculations.

The Poisson-Boltzmann screening is used, because in that case, the potential at far away place is defined as zero, so the absolute Fermi energy is well defined. This option should be used with care. The E\_Fermi is the Fermi energy in the unit of eV (usually is negative). Note, in this case, the total number of electron NUM\_ELECTRON will be adjusted automatically, and its input value will only be used as an initial value. mix\_Q is a charge mixing parameter for the total charge. Small value will be more stable, but slower. Suggest to try 0.1. The drho\_pulay is the value for charge density selfconsistent error drho before to turn on the pulay charge mixing. The pulay should not be used at the beginnning, that will cause instability. Suggest to try drho\_pulay at 0.06. Unfortunately, one might has to adjust mix\_Q, drho\_pulay for better convergence and speed. Also, one might increase the smearing of Fermi-Dirac occupation function (in the line of SCF\_ITER0) to make the calculation stable (e.g., 0.25 eV). Besides, one usually only do a JOB=SCF calculation, with FIX\_FERMI=T, instead of doing RELAX. For relaxation, one can first do a SCF with FIX\_FERMI=T, then do RELAX with FIX\_FERMI=F, but use the total charge output (reported in REPORT) as NUM\_ELECTRON to do a RELAX calculation. One can iterate this loop.

### CONSTRAINT\_MAG

|Tag|CONSTRAINT\_MAG|
| --- | --- |
|**Format**|CONSTRAINT\_MAG = 0 / 1|
|**Default**|CONSTRAINT\_MAG = 0|

This input line is used to put constraint on magnetic moment at each atom for spin=2 calculations. If CONSTRAINT\_MAG=1, this will be turned on, if CONSTRAINT\_MAG=0, this will not be used. Default is CONSTRAINT\_MAG=0. Currently, it only works for spin=2. This is to add an energy penalty term:

$$\sum_{iatom} alpha\_mag(iatom) (M(iatom)-M_{input}(iatom))^2$$

in the total energy expression. Note, M(iatom) is the magnetic moment calculated from spin up and down charge density using Hirshfield method, much like in the CHARGE\_DECOMP. M\_{input}(iatom) are the input magnetic moment in atom.config file, under a section name: CONSTRAINT\_MAG. alpha\_mag(iatom) is also input from the atom.config from that CONSTRAINT\_MAG section. Note, if this section is not provided in atom.config, the default value for M\_{input}(iatom) is zero, and alpha\_mag(iatom) is also zero. The unit of alpha\_mag is in eV.

Alpha\_mag should not be set to be too large. Otherwise the SCF iteration will be difficult to converge. One recommended alpha\_mag is 0.01 eV. One can also set different alpha\_mag for different atoms. Also, note that one can use the above additional energy to add an effective magnetic field H at each atom (with different amplitude etc). To do that, one can use a large M\_{input}(iatom) in the direction one wants, in the mean time, uses a very small alpha\_mag(iatom) which is proportional to 1/M\_{input}(iatom). That will provide an effective H field on each atom.

### SPIN222\_MAGDIR\_STEPFIX

|Tag|SPIN222\_MAGDIR\_STEPFIX|
| --- | --- |
|**Format**|SPIN222\_MAGDIR\_STEPFIX = N|
|**Default**|SPIN222\_MAGDIR\_STEPFIX = 0 (for XCFUNCTIONAL = LDA)|
||SPIN222\_MAGDIR\_STEPFIX = 1000 (for everything else)|

This input line is used to fix the direction of magnetic moment after SPIN222\_MAGDIR\_STEPFIX self-consistent iterations. If SPIN222\_MAGDIR\_STEPFIX = 30, it means that after 30 self-consistent interations, the direction of magnetic moment will be fixed. But for XCFUNCTIONAL = LDA, default SPIN222\_MAGDIR\_STEPFIX is 0, the direction of magnetic moment will not be fixed. If you set SPIN222\_MAGDIR\_STEPFIX = 1, initial direction of magnetic moment will be fixed, and will not change with self-consistent interations. Please note that SPIN222\_MAGDIR\_STEPFIX is only used for SPIN = 222.

### E\_FINITE

|Tag|E\_FINITE|
| --- | --- |
|**Format**|E\_FINITE = T / F  Ex Ey Ez|
|**Default**|E\_FINITE = F 0 0 0|

This parameter is used to set the homogeneous electric field in the electric enthalpy functional. If the first parameter is T, it will compute the self-consistent response to finite electric fields. Ex,Ey,Ez in unit eV/Angstrom.

In the output of screen, one can search for "Pel", the electronic dipole moment, in unit e*Angstrom.The first "Pel" is the result with fields turned off. The second "Pel" is the result with fields turned on.

An example of AlAs, the file atom.config:

```bash
 2
 LATTICE
 4.0543775558         0.0000000000         0.0000000000
 2.0271887779         3.5111939599         0.0000000000
 2.0271887779         1.1703979866         3.3103854121
 POSITION
 13 0.000000000       0.000000000         0.000000000  1 1 1
 33 0.749999998       0.750000011         0.749999994  1 1 1
```

the file etot.input:

```bash
4  1
JOB = scf
IN.PSP1 = Al.SG15.PBE.UPF
IN.PSP2 = As.SG15.PBE.UPF
IN.ATOM = atom.config

Ecut    = 50
Ecut2   = 200
MP_N123 = 4 4 4  0 0 0 2

e_finite  = T 0.00 0.00 0.001

precision = double
e_error   = 0.0
wg_error  = 0.0
rho_error = 1.d-6

out.force = t
 
```

The OUT.FORCE=T is needed to calculate born effective charge, for accurate forces, Ecut2=4*Ecut is recommended. Kpoints are set by MP\_N123 without symmetry(E\_FINITE=T can not use symmetry). Kpoints grid should be large enough to get converged results. 
    
PRECISION=DOUBLE is set for more accurate results, much more here we set e\_error=0, wg\_error=0, rho\_error=1.d-6 to ensure the wavefunctions in good convergency.

For systems with small gap, one need to use smaller FremidE, also need to use much smaller strength of field in case of no local minimum, like following example, the file atom.config:
      
```bash
2
LATTICE
     4.06599283     0.00000000     0.00000000
     2.03299642     3.52125308     0.00000000
     2.03299642     1.17375103     3.31986925
POSITION
 31     0.00000000     0.00000000     0.00000000 1 1 1
 33     0.25000000     0.24999999     0.25000000 1 1 1

```

the file etot.input:

```bash
4  1
JOB = scf
IN.PSP1 = Ga.SG15.PBE.UPF
IN.PSP2 = As.SG15.PBE.UPF
IN.ATOM = atom.config

Ecut    = 50
Ecut2   = 200
# for correct and converged results, maybe need to use more kpoints
MP_N123 = 4 4 4  0 0 0 2   
e_finite  = T 0.00 0.00 0.0001   

precision = double
e_error   = 0.0
wg_error  = 0.0
rho_error = 1.d-6

out.force = t

#use FermidE=0.001eV
SCF_ITER0_1 =    6   4    3    0.0000     0.0010    1   
SCF_ITER0_2 =   94   4    3    1.0000     0.0010    1 

```

### RVV10\_DETAIL

|Tag|RVV10\_DETAIL|
| --- | --- |
|**Format**|RVV10\_DETAIL = b, c|
|**Default**|RVV10\_DETAIL = 6.3, 0.0093|

The current code implemented the RVV10 dispersion interaction (G. Roman-Perez, J.M. Soler, Phys. Rev. Lett. 103, 096102 (2009); R. Sabatini, T. Gorni, S. de Gironcoli, Phys. Rev. B, 87, 041108 (2013)). The RVV10 use a nonlocal integral of charge densities at different points, r, r' derived from RPA formalism. The kernel of this integral is  

$$\phi^{VV10}(r,r')=-{3e^4\over 2m^2} {1\over g g' (g+g')}$$

Here $g=\omega_0(r)(r-r')^2+k(r)$, $g'=\omega_0(r')(r-r')^2+k(r')$. Furthermore, $\omega_0(r)=\sqrt{\omega_g^2(r)+\omega_p^2(r)/3}$. $\omega_p^2(r)=4\pi{n(r)}e^2/m$ is the plasma frequency.

$\omega_g^2(r)=c(\hbar^2/m^2)|\frac{\nabla{n(r)}}{n(r)}|^4$ and $k(r)=3 \pi b({n(r)\over9\pi})^{1\over 6}$, Here b and c are parameters. For default, b=6.3, c=0.0093. But one can use RVV10\_DETAIL to specify different b and c values. This is only used when XCFUNCTIONAL contains RVV10.

## I/O tags

### IN.ATOM

|Tag|IN.ATOM|
| --- | --- |
|**Format**|IN.ATOM = atom.config|
|**Default**|NO DEFAULT|

IN.ATOM is used to read the atomic positions file, this file contains the lattice geometry and ionic positions, optional tags -- force, velocity, magnetic, constraint\_mag, magnetic\_xyz, langevin\_atomfact\_tg, stress\_mask, et al. Its specification is described in the section \ref{inputfile:atomconfig} of this manual.

### IN.PSP

|Tag|IN.PSP|
| --- | --- |
|**Format**|IN.PSP1 = H.NCPP.UPF|
||IN.PSP2 = C.NCPP.UPF|
||...|
|**Default**|NO DEFAULT|

The names of the pseudopotential files. `IN.PSP1' is the first atom type, `IN.PSP2' is the second atom type, and the rest can be deduced by analogy. The order of different element types is arbitrary. Please see section \ref{inputfile:pseudopotential} for a discussion of different pseudopotentials.

### IN.KPT

|Tag|IN.KPT|
| --- | --- |
|**Format**|IN.KPT = T / F|
|**Default**|IN.KPT = F|

IN.KPT = T, PWmat will use the k-points from file `IN.KPT` which contains the k-points and their weights. %The IN.KPT can be generated (together with IN.SYMM) by running `check.x` with information from variable MP\_N123. Note, IN.KPT, IN.SYMM usually work together. IN.KPT has a higher priority than MP\_N123. If IN.KPT = F, PWmat will not use the file `IN.KPT`. PWmat will always output koints in `OUT.KPT` file. Please check the IN.KPT(OUT.KPT) subsession for more detail format about this file.

### IN.SYMM

|Tag|IN.SYMM|
| --- | --- |
|**Format**|IN.SYMM = T / F|
|**Default**|IN.SYMM = F|

IN.SYMM = T, PWmat will use the file `IN.SYMM` (the name is fixed) to perform symmetry operations. The PWmat supports space group symmetry for crystals. `IN.SYMM` should contain space group symmetry operations. %`IN.SYMM` is usually generated (together with IN.KPT) by running `check.x` (which will also check whether the IN.SYMM exists if IN.SYMM=T). Usually, symmetry can be generated automatically by using MP\_N123 line. However, one can also copy over the previous OUT.SYMM into IN.SYMM for explicit input (e.g., one can even delete some symmetry operations). The IN.SYMM usually should work together with IN.KPT (for the reduced k-points). Note, if both IN.SYMM=T, IN.KPT=F are specified, and also MP\_N123 are also specified, PWmat will generate Kpoints using the symmetry operations provided in file `IN.SYMM`. If IN.SYMM=T, IN.KPT=T and also MP\_N123 are also specified, the IN.KPT=T has a higher priority, MP\_N123 will not be used. If IN.SYMM = F, PWmat will not use file `IN.SYMM`, This is the default value. PWmat will always output symmetry operations in `OUT.SYMM` file. Please check the IN.SYMM(OUT.SYMM) subsession for more detail format about this file.


### IN.OCC

|Tag|IN.OCC|
| --- | --- |
|**Format**|IN.OCC = T / F|
|**Default**|IN.OCC = F|

Related items: PROJ3\_DETAIL, IN.iproj3\_CC\_2spin, IN.OCC\_T, IN.CC. In all these options, the PWmat will try to use special ways to determine the occupation of wave functions, or eigen states, instead of using the conventional Fermi-Dirac distribution from the SCF\_ITER0\_1, SCF\_ITER0\_2, SCF\_ITER1\_1 lines. Note, in general, not necessarily the eigen states will be occupied, instead an linear combination of the eigen states can be occupied. These options will be particularly useful for either constraint DFT (with some excited  electron states, or empty hole states) SCF or RELAX jobs, or TDDFT simulations. For example, it can be used to prepare some special initial excited state in TDDFT.

In this option, PWmat will read a file called `IN.OCC`. This is to specify the occupation for each Kohn-Sham orbital $\phi_j$ for a constraint DFT calculation. 

>
>**WARNING**: `IN.OCC=T` is equivalent to `IN.OCC=T,0`
{: .block-warning}

In the following, we will use $\phi_i$ to denote the adiabatic eigen state of the Kohn-Sham Equation: $H \phi_i = \epsilon_i \phi_i$. In the meantime, we can input another set of wave function: $\{\psi_j\}$, e.g, through the IN.WG=T option. This $\{\psi_j\}$ is a bit like the time evolving wave functions in TDDFT, and they can be used to help the
occupation of states to generate the charge density.

The following are the different options for iproj.

**iproj=0**: (IN.OCC = T, or IN.OCC = T, 0): The SCF calculation charge density will be generated as:  $\rho(r)= \sum_i o(i) |\phi_i(r)|^2$, and the occupation number $o(i)$ will be input from the IN.OCC file (see the explanation below). Note, for iproj=0, the o(i) is fixed (for which state is which). For example, if the third state in IN.OCC is empty (o(3)=0), then during the SCF calculation, or atomic relaxation, it is always the third adiabatic state which is unoccupied. This is okay for simple cases, but for more complicated cases, it can cause problem, since during the SCF iteration, or atomic relaxation, the order (which state is the third) according to the adiabatic state eigen energies can often change (re-ordered). As a result, the third state might not be the physical state you like to keep it empty. One option in that case is to use iproj=1.

**iproj=1**: In this option, the index of which state is which is determined by a projection between the adiabatic eigen state $\phi_i$ and the input state $\psi_j$ from IN.WG. So, in order to use this, one has to have IN.WG. More specifically, map(i) equals the j which maximizes $|<\phi_i|\psi_j>|^2$. Basically, this identifies which $\phi_i$ is the input $\psi_{map(i)}$. As a result, the charge density is generated as: $\rho(r)= \sum_i o(map(i)) |\phi_i(r)|^2$. Note, { $\psi_j$ } is input from IN.WG. Once again, o(j) is input from IN.OCC, which has the following form (for iproj=0,1,3):

```bash
o1,o2,o3,.....o_mx        ! for kpt=1
o1,o2,o3,.....o_mx        ! for kpt=2
........
o1,o2,o3,.....o_mx        ! for kpt=nkpt
```

Here mx is the num\_band, and nkpt is the number of reduced kpoint. Thus, in each line, there are mx number, and there are nkpt line. Note, one can write something like: 4*1,2*0 to replace: 1,1,1,1,0,0. Note, o(j,kpt) is the occupation number discussed above (it must be between 0 and 1, even for SPIN=1, should not be 2). Note, the full occupation is o(j,kpt)=1. So, if SPIN=1, o(j,kpt)=1 means this orbital will occupy 2 electron, and o(j,kpt)=0.5 means this orbital will occupy 1 electron.

If SPIN=2, then one has also to provide an file: IN.OCC\_2, which has the same format, but specify the spin down component occupation. In that case, o(j,kpt)=1 means this orbital (up or down spin) will occupy one electron, and o(j,kpt)=0 means
this orbital of this spin will occupy zero electron.

**iproj=2 (or 22)**:  in the above case of iproj=0 and 1, the Fermi-Dirac distribution function will be ignored (not used, or it is like Fermi-Dirac=0). But even for iproj=1, there could be cases where it is difficult to identify which $\phi_i$ is $\psi_j$, for example, if several states have the similar amplitude overlaps (e.g., around 0.3-0.5). In that case, even if we select and occupy one $\phi_i$, the results might not be good either. What we really need is to construct one $\psi_j$-like wave function $\psi'_j$ from a linear combination of $\phi_i$, then occupy $\psi'_j$. More specifically, we can have: $\psi'_j = \sum_i f(i,j) <\phi_i|\psi_j> \phi_i$. Here f(i,j) is a selection and orthonormalization factor. For iproj=22, there is no selection, and f(i,j) is just for orthonormalization (among all the $\{ \psi'_j \}$). For iproj=2, a selection weight factor is used for different i. More specifically, $f(i,j)=exp(-((1-x)/0.7)^6)$, where $x=|<\phi_i|\psi_j>|^2$. Thus, effectively, it only select the $\phi_i$ states with overlap larger than 0.3. After this selection factor, it is orthonormalized. The idea here is only to select a few $\phi_i$, and use their linear combination to construct a state resemble that of the original $\psi_j$ state, but not to completely reconstruct the $\psi_j$ state. This will be useful to deal with the case where two states anticross each other, so the $\phi_i$ character has changed, but a linear combination can reconstruct the original $\psi_j$.

In terms of occupation, it uses a different strategy than iproj=0,1,3. For iproj=2, the Fermi-Dirac distribution is still used, then on top of it, an exception for some band is used to add or substract some states. The idea is to use this to add one excited electron or subtract one hole  (on $\psi_j$) (called exception states) from the otherwise Fermi-Dirac calculation. Note, when decide the Fermi energy, a total charge of NUM\_ELECTRON - dcharge is used, and dcharge is the charge from $\psi_j$ as specified in the IN.OCC to be explained below. So, the actual total charge is NUM\_ELECTRON. Also note, not only the charge density of $\psi_j$ is contributed, also its kinetic and nonlocal potential part of the energy. The total charge density can be specified as: $\rho(r)=\rho_{FD}(r)+ \sum_j o(j) |\psi'_j(r)|^2$, here $\rho_{FD}(r)$ is the Fermi-Dirac charge density, and o(j) is the occupation from the IN.OCC, which has  the following form in the case of iproj=2(or 22):

```bash
o1,o2,o3,.....o_mx        ! for kpt=1
o1,o2,o3,.....o_mx        ! for kpt=2
........
o1,o2,o3,.....o_mx        ! for kpt=nkpt
-------- (above nkpt lines are not used, this dashed line must be here)
nump                          ! The number of exception states
iband(1,1),od(1,1)            ! for indx=1,   kpt=1
iband(2,1),od(2,1)            ! for indx=2,   kpt=1
......
iband(nump,1),od(nump,1)      ! for indx=nump,kpt=1
iband(1,2),od(1,2)            ! for indx=1,   kpt=2
iband(2,2),od(2,2)            ! for indx=2,   kpt=2
......
iband(nump,2),od(nump,2)      ! for indx=nump,kpt=2
......                        ! repeat for different kpoint
......
iband(1,nkpt),od(1,nkpt)            ! for indx=1,   kpt=nkpt
iband(2,nkpt),od(2,nkpt)            ! for indx=2,   kpt=nkpt
......
iband(nump,nkpt),od(nump,nkpt)      ! for indx=nump,kpt=nkpt

```


Note, in this file, there must be nump\*nkpt lines after the first nump+1 lines. The nump is the number of exception states in each kpoints. Each of the nump\*nkpt line should have two numbers. The first number is the band index, the second number is an occupation (can be positive, or negative (hole)).

Note, for SPIN=2, one needs to have the same format for IN.OCC\_2.

Note, iproj=0,1,2 are often used for both SCF calculation, as well as RELAX. iproj=0 is also often used
for TDDFT calculation.

For the cases of iproj=1 for  RELAX, the $\psi_{j}$ for next RELAX step is updated from the $\phi_{map(j)}$ from the previous relaxation step. For the case of iproj=2 for RELAX, $\psi_j$ for next RELAX step is replaced with $\psi'_j$ from previous RELAX step. Doing this will allow one to relax the system even if there is a state crossing for a hole state (or electron state), while tracking and keeping the same hole or excited electron.

The iproj=2 is usually more stable than iproj=1 for both SCF and RELAX. So, it should be used if iproj=0 and 1 are unstable, and do not converge.

**iproj=3 (or 33)**:  this is used for SCF or TDDFT calculations, not for RELAX. In iproj=2, the exception is used to calculate the charge density. It is mostly designed to deal with RELAX, and conceptually, we like to eventually make the $\psi_j$ the adiabatic eigen states. There are cases where we are interested in occupy specific input wave function $\psi_j$, but they are not the adiabatic eigen states $\phi_i$. More importantly, the occupied state $\psi'_j$ should be completely represented by $\phi_i$. So, the occupied state is not the exact $\psi_j$ input from IN.WG. If the fixed $\psi_j$ is to be occupied, one can use JOB=WKM. Here, we will occupy $\psi'_j = \sum_i f(i) <\phi_i|\psi_j> \phi_i$, then orthonormalized following this formula. Thus, $\rho(r)=\sum_j o(j) |\psi'_j(r)|^2$,  here o(j) is input from IN.OCC. In the formula, f(i) is used to provide some possible truncation (so there is no high energy $\phi_i$ components for $\psi'_j$. This can be provided by the following input line:

```bash
PROJ3\_DETAIL = Ecut\_proj3, dEcut\_proj3  (in unit of eV)
```

Then: $f(i)=1/(exp((\epsilon(i)-Ecut\_proj3)/dEcut\_proj3)+1)$

If the PROJ3\_DETAIL does not exist, then f(i)=1

Also, in the IN.OCC, if one o(j) is negative, then $\psi'_j=\phi_j$, and $o(j)=|o(j)|$.

This can be used to provide the initial state for TDDFT calculation (can be used together with JOB=TDDFT). So, the wanted initial TDDFT $\psi_j(t=0)$ wave function (whatever wave function, not necessarily eigen states) can be input from IN.WG. But for our TDDFT implementation, we need to gaurantee the $\psi_j$ can be represented by the set of eigen states $\{ \phi_i \}$. This IN.OCC and the first step in TDDFT (or the SCF calculation) can just gaurantee this through a self-consistent iteration. Note, $\psi_j$ can be some localized state, or some ionic states for the ion far away from a colliding surface. One can use other means to pre-construct IN.WG. Note, this iproj=3 will automatically output a file OUT.iproj3\_CC\_2spin, it contains the coefficient: $CC(i,j)=<\phi_i|\psi_j>$ with orthonormalization. For iproj=33, this OUT.iproj3\_CC\_2spin must be copied into IN.iproj3\_CC\_2spin, and inside etot.input, one has to place one line: IN.iproj3\_CC\_2spin=T. In this case, the IN.WG input is actually the eigen state $\phi_i$ instead of $\psi_j$, and $\psi_j$ is constructed as: $\psi_j=\sum_i CC(i,j) \phi_i$. In this way, one can input both the eigen state and the $\psi_j$, thus continue a iproj=3 calculation.

**Comparison with other options**: Note, for TDDFT calculation, one can compare IN.OCC=T,3 calculation with IN.OCC\_T=T calculation. In IN.OCC\_T, the occupation of each state o(j) can be changed with time, specified inside IN.OCC\_T. That is a powerful tool to actually change the occupation (e.g., to describe one electron gradually disappear due to some other physical process). It can also be used to prepare the initial state in a TDDFT calculation, e.g, quickly remove one state. But this is less general than IN.OCC=T,3, it can also has some stability issues. Also note, if one very quickly remove one electron through IN.OCC\_T, the physical meaning might be different from the initial condition propered using IN.OCC=T,0. In IN.OCC=T,0, the initial condition is in a equilibirium condition, but that is not the case for IN.OCC\_T calculation.

Lastly, one can also compare IN.OCC=T,3 with IN.CC. IN.CC can also be used to prepare a mixing state as the initial state for TDDFT calculation. It is easy to construct the initial wave function, but it is less general, and less powerful than the IN.OCC=T,3 procedure.

### IN.OCC\_T

|Tag|IN.OCC\_T|
| --- | --- |
|**Format**|IN.OCC\_T = T|
|**Default**|IN.OCC\_T = F|

If IN.OCC_T = T, PWmat will read in file `IN.OCC\_T'  (and `IN.OCC\_T\_2' if spin=2), and it will ignore the imth\_Fermi  flag In the SCF\_ITERx\_x line.

### IN.NONSCF

|Tag|IN.NONSCF|
| --- | --- |
|**Format**|IN.NONSCF = T/F|
|**Default**|IN.NONSCF = F|

This parameter is used to set optional NONSCF parameter in file IN.NONSCF.

### IN.RELAXOPT

|Tag|IN.RELAXOPT|
| --- | --- |
|**Format**|IN.RELAXOPT = T/F|
|**Default**|IN.RELAXOPT = F|

This parameter is used to set optional RELAX parameter in file IN.RELAXOPT.

### IN.MDOPT

|Tag|IN.MDOPT|
| --- | --- |
|**Format**|IN.MDOPT = T/F|
|**Default**|IN.MDOPT = F|

This parameter is used to set optional MD parameter in file IN.MDOPT.
    
If the method of MD is 2,3,4 or 5, one can use the file IN.MDOPT to set detailed parameters by setting IN.MDOPT=T. See more in \ref{otherinput:in.mdopt}.

### IN.TDDFTOPT

|Tag|IN.TDDFTOPT|
| --- | --- |
|**Format**|IN.TDDFTOPT = T/F|
|**Default**|IN.TDDFTOPT = F|

This parameter is used to set optional TDDFT parameter in file IN.TDDFTOPT.

### IN.EXT\_FORCE

|Tag|IN.EXT\_FORCE|
| --- | --- |
|**Format**|IN.EXT\_FORCE= T / F|
|**Default**|IN.EXT\_FORCE = F|

This parameter is used to provide an external force (unit eV/amstrong) for each atom during MD simulation. See more in MD\_DETAIL. If IN.EXT\_FORCE=T, a file IN.EXT\_FORCE will be provided, it has the following format:

```bash
    natom
    iatom, fx, fy, fz    ! unit eV/Amstrong
    .................
    iatom, fx, fy, fz    ! There will be natom lines
```

### IN.SOLVENT

|Tag|IN.SOLVENT|
| --- | --- |
|**Format**|IN.SOLVENT = T / F|
|**Default**|IN.SOLVENT = F|

When calculating the energy of a solute molecule in a liquid solvent, e.g., in electric chemistry study, there are two possible approaches. One is to use explicit solvent molecule (e.g., water molecules) and carry out molecular dynamics simulations, another is to use implicit solvent models. The implicit solvent model represents the effect of the solvent with a continuum mediate, mostly includes its effects of electric static polarization. Compared with the explicit solvent molecules and molecular dynamics, the implicit solvent model is much faster. We have followed the work of self-consistent continuum solvation (SCCS) model \cite{solvent1}, as well as similar formalism in Ref.\cite{solvent2}. In quantum chemistry, this is also called: polarizable continuum model (PCM). We have also implemented the linearized Poisson-Boltzmann screening for the effects of free ions (salt, or H+, OH- in low or high pH value situations) \cite{solvent3}. These models use a continuum mediate and a space variation dielectric constanr $\epsilon(r)$ to represent the solvation effects, mostly the polarization effects. There are three energy terms: the polarization solvation energy, the cavity energy (the surface tension, or can also be considered as surface van der Waals energy), and volume energy (pressure energy, PV). If Poisson-Boltzmann equation is used to describe the ion screening, there is another term which describes the ion energy within an electric static potential. When the Poisson-Boltzmann equation is used, the absoulte potential zero is defined at the far away place (there is no ambiguity of the absolute potential position). This also allows us to define the absolute potential of the electrode (for example, the standard hydrogen electrod, SHE, potential in water is at -4.42 eV). In this situation, we can use a fixed potential (fixed Fermi energy calculation). This is controlled using Fix\_Fermi = T in etot.input (please see the corresponding item in the manual).

**IN.SOLVENT = T**, will use solvent model. Note, all the other input in etot.input will be the same. In other words, solvent model can be used to do single SCF, RELAX, MD and TDDFT calculations. However, one has to prepare a `IN.SOLVENT` file in the running directory, which control the parameters for the solvent model. See more in \ref{otherinput:in.solvent}.

When IN.SOLVENT=T, after the PWmat run, a OUT.SOLVENT file will be generated which lists all the options used for the solvent model.


### IN.A\_FIELD

|Tag|IN.A\_FIELD|
| --- | --- |
|**Format**|IN.A\_FIELD= T / F, a\_field1, a\_field2, a\_field3|
|**Default**|IN.A\_FIELD= F 0.0 0.0 0.0|

>
>**TIP**: for PWmat version later than 20200824, it has a new format:
>
> IN.A\_FIELD\_LIST1= a\_field1, a\_field2, a\_field3 IN.TDDFT\_TIME1
> IN.A\_FIELD\_LIST2= a\_field1, a\_field2, a\_field3 IN.TDDFT\_TIME2
> ... 
> the maximum support is 20 rows. This can be used to add a circularly polarized light.
> IN.TDDFT\_TIME1,IN.TDDFT\_TIME2... is the name of TDDFT\_TIME file. You need to prepare same number of TDDFT\_TIME files, it has the same format as IN.TDDFT\_TIME,
>
>        0 ftddft(0)
>        1 ftddft(1)
>        ...
>        N ftddft(N)

This controls the G-space external potential input for tddft calculation(only used when TDDFT\_SPACE=-1,...).The tddft hamiltonian,

$$
H=1/2 (-i\nabla_x + a\_field1)^2+1/2(-i\nabla_y+ a\_field2)^2+1/2(-i\nabla_z + a\_field3)^2
$$

 The values of $a\_field1,2,3$ are all in units of 1/Bohr.

 ### IN.WG


|Tag|IN.WG|
| --- | --- |
|**Format**|IN.WG = T / F|
|**Default**|IN.WG = F|

IN.WG = T, PWmat will read in the initial wave functions in G-space from the file `IN.WG` (e.g., from previous calculation, copied over from OUT.WG, but note the node1 from the current calculation must be the same as in the previous calculation to generate OUT.WG). When SPIN = 2, an extra file `IN.WG\_2` will also be read in. Note, IN.WG, OUT.WG can be plotted using utility program `plot\_wg.x`, which can be used to view each wave function.  For people like to see the format of IN.WG, OUT.WG, one can check the plot\_wg.f90 file which is source code of `plot\_wg.x` utility program. If IN.WG = F, the PWmat will start from random wave function.

### IN.RHO

|Tag|IN.RHO|
| --- | --- |
|**Format**|IN.RHO = T / F|
|**Default**|IN.RHO = F|

IN.RHO = T, PWmat will read in the initial charge density from file `IN.RHO`, stored in the real space grid (N1L, N2L, N3L). This can be copied over from OUT.RHO of previous calculation. Note, the node1 in current calculation, and previous calculation to generation OUT.RHO must dividable from one way or the other. Note, if both IN.VR and IN.RHO are set to T, the program will use the read-in potential to start the calculation. If SPIN = 2, PWmat will read an extra file `IN.RHO\_2`. IN.RHO=T is also needed for JOB=POTENTIAL. If SPIN = 22, only a single IN.RHO will be needed. If SPIN = 222, besides IN.RHO, a IN.RHO\_SOM (a complex 2x2 spin matrix density) will be needed. If IN.RHO = F, not input the charge density.

    One can use utility program: convert\_rho.x to plot OUT.RHO or IN.RHO. One can also check convert\_rho.f90 for the format of IN.RHO, OUT.RHO.

### IN.RHO\_ADD

|Tag|IN.RHO\_ADD|
| --- | --- |
|**Format**|IN.RHO\_ADD = T / F|
|**Default**|IN.RHO\_ADD = F|

IN.RHO\_ADD = T, PWmat will read in an additional charge density $\rho_{add}$ from input file: IN.RHO\_ADD, and this $\rho_{add}$ will be added to the total charge density calculated from the wave functions. In another word, $\rho(i)=\sum_i |\psi_i(r)|^2 occ(i) + \rho_{add}(r)$. Note, this $\rho_{add}$ will not be counted as part of NUM\_ELECTRON when determing $occ(i)$. This function can be used for many different algorithms. In particular, JOB=SCFEP is one particular case of this (but please continue to use JOB=SCFEP). Note, one can use the utility file: convert\_wg2rho.f to general IN.RHO\_ADD file from the output wave function file OUT.WG. When spin=2, one also needs to provide an IN.RHO\_ADD\_2 file for spin down additional charge density. If IN.RHO\_ADD = F, no additional charge density is used.

### IN.VR

|Tag|IN.VR|
| --- | --- |
|**Format**|IN.VR = T / F|
|**Default**|IN.VR = F|

IN.VR = T, PWmat will read in the initial potential from file `IN.VR` in real space grid: (N1L, N2L, N3L). Note, if both IN.VR and IN.RHO are set to T, the program will use the read-in potential to start the calculation. The format, and requirement of IN.VR are the same as that for IN.RHO. When SPIN = 2, an extra file `IN.VR\_2` will also be read. When SPIN=22, just a single IN.VR will be read (no IN.VR\_2). When SPIN=222, besides IN.VR, IN.VR\_SOM (complex 2x2 spin matrix potential), and IN.VR\_DELTA (a single real up-down potential) need to be read in. If IN.VR = F, not read in the file.

### IN.VEXT

|Tag|IN.VEXT|
| --- | --- |
|**Format**|IN.VEXT = T / F|
|**Default**|IN.VEXT = F|

IN.VEXT = T, PWmat will read in an external potential from file `IN.VEXT` in real space grid (N1L, N2L, N3L). Both the total energy and forces are calculated using this external potential. This can be useful to calculate the influence of an external potential (e.g., an electric field), for JOB=SCF, RELAX or MD. Note, the IN.VEXT has the same format as that in IN.RHO and IN.VR. Its unit is Hartree. One can try (and check, modify) the utility programs: calculate\_Vext.f90 and gen\_external\_efield.f90 to generate such IN.VEXT. One can also write such IN.VEXT by self-made codes. If IN.VEXT = F, no external potential is used.

### IN.LDAU

|Tag|IN.LDAU|
| --- | --- |
|**Format**|IN.LDAU = T / F|
|**Default**|IN.LDAU = F|

IN.LDAU = T, PWmat will read in the initialization of LDA+U from file `IN.LDAU` (Note this is not the U parameters, instead they are the LDA+U occupation values calculated from JOB=SCF calculation). This setting is only required when JOB=NONSCF.  One needs to copy OUT.LDAU to IN.LDAU after the JOB=SCF. If SPIN = 2, PWmat will read an extra file `IN.LDAU\_2`, also there will be OUT.LDAU\_2 after JOB=SCF. If IN.LDAU = F, not input the initialization of LDA+U.

### OUT.WG

|Tag|OUT.WG|
| --- | --- |
|**Format**|OUT.WG = T / F|
|**Default**|OUT.WG = T|

OUT.WG = T, PWmat will output a file `OUT.WG`, which stores the final wave functions in G-space. When SPIN = 2, an extra file `OUT.WG\_2` will also be output. This is the default value. Use utility program  `plot\_wg.x` to plot the wave functions. More details about `plot\_wg.x`, please refer to PWmat website \url{http://www.pwmat.com/utility-download}.

If OUT.WG = F, will not output the wave function file.

### OUT.RHO

|Tag|OUT.RHO|
| --- | --- |
|**Format**|OUT.RHO = T / F|
|**Default**|OUT.RHO = T|

OUT.RHO = T, PWmat will output a file `OUT.RHO`, the final charge density in real space grid (N1L, N2L, N3L). This is the default value. If SPIN = 2, PWmat will write out an extra file `OUT.RHO\_2`. If SPIN=222, PWmat will also output OUT.RHO\_SOM, a 2x2 complex spin matrix density.

Use utility program `convert\_rho.x` to plot the charge density (unit in e/Bohr$^3$). More details about `convert\_rho.x`, please refer to PWmat website http://www.pwmat.com/utility-download.

If OUT.RHO = F, not output the charge file.
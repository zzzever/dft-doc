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
    
If **iflag\_wkm\_Hxc=1**, then the bands: $[1, nb\_exclude\_Hxc]$ counted from the bottom will not be included in the exchange-correlation function evaluations, and in this case the **ss1\_u** and **ss2\_u** are used to occupy the Wannier function so to get a closed shell structure.

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

>WARNING
>
>It is notallowed to use this parameter when ``ENERGY\_DECOMP = T'' in etot.input. If NUM\_BLOCKED\_PSI = T,the decomposed energy can suddenly be very wrong.
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

$$$H_{m,m'}(T) = \int dr \chi^*_m(r-\tau_m) \hat{H} \chi_{m'}(r-(\tau_{m'}+T))$$, 

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

### ECUT2

### ECUT2L

### ECUTP

### N123
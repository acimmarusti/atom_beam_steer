
# Atomic beam steering simulation

This ROOT macro aims to simulate and analyze the propagation of an atomic beam through lasers (forming a 2-D optical molasses) and into a designated target. It has been developed since spring 2010 for the purpose of simulating a proposed upgrade to our cavity quantum electrodynamics experiment.

## NEWS:

[2013-03-15] The code has been tested on ROOT versions 5.26 through 5.34 only on GNU/Linux workstations (CentOS and Debian).

[2011-07-01] Our paper using this software has been published in the Revista Mexican de Física, as part of the Quantum Optics conference proceedings: http://rmf.smf.mx/pdf/rmf-s/57/3/57_3_29.pdf. For more info see publication list.

## AUTHORS:
The macros were developed by Andres D. Cimmarusti. The code has benefited from discussions with David G. Norris and Luis A. Orozco.

## DESCRIPTION:

The simulation starts by creating an atomic beam from known experimental results. We chose a Low Velocity Intense Source (LVIS) of atoms as the example. The paper by Lu et al 10.1103/PhysRevLett.77.3331, characterizes a type of atomic beam.

In an LVIS, the atoms escape a Magneto-optical trap through a hole in a mirror. We initially generate the random (uniformly distributed) (x,y) position of an atom as it exits such a hole. Using a geometric argument involving the distance between and hole and the center of the atomic trap, we randomly assign the polar angle of the atom's trajectory. Its velocity is again randomly generated from a Gaussian distribution as explained in the paper cited above. This concludes the generation of our initial atomic beam.

The atoms propagate freely under the influence of gravity until they reach a 2D optical molasses. The optical molasses consists of a three dimensional region of space where four Gaussian laser beams intersect. This was possible thanks to the heavy use of ROOT's geometry classes. The geometry, including the Gaussian laser beam parameters, admit extensive tweaking. This part of the program contains a Monte Carlo method which we use to simulate the interaction of the lasers with each atom (photon scattering and subsequent spontaneous emission). The overall effect of the 2D optical molasses on the atoms is to provide a damping force to slow down their average speed in the plane of the molasses. This is, by far, the most resource intensive chunk of the code. Despite that, the macro (if compiled) can run on modest hardware (Intel Core 2 Duo 2.66 GHz) for 100000 atoms in about 10 minutes.

The atomic beam is effectively "bent" or "steered" upong leaving the 2D optical molasses region. It continues to propagate freely until it reaches its goal. The goal in the original simulation of our experiment was the spatial TEM00 mode formed in a high-finesse optical cavity. In practice, this mode was about 56 µm in diameter and with a length of the separation between the mirrors of the cavity. In this phase, it can estimate the fraction of atoms in the beam that eventually passed through the cavity mode region.

Alternatively after the 2D optical molasses, the macro has the capabilities of implementing a so-called laser lens. The idea was to use highly focused, counterpropagating Gaussian laser beams to induce 2D spatial confinement of the atomic beam and into the cavity mode. Enabling this part of the code causes a significant performance penalty, as the code is esentially similar to the 2D optical molasses, except that the beams are highly focused to provide the necessary position dependent force on the atoms via the intensity gradient. Using the extra feature can boost the fraction of atoms entering the cavity mode by a factor of 2.

## LICENSE:

This program is released under the terms of the GNU Public License v3 (GPLv3). For more information please see:

http://www.gnu.org/licenses/gpl-3.0.html

## HOW TO USE:

Before trying to run the script, it is necessary to have ROOT installed on your system. ROOT is a C++ framework. As such, it brings lots of scientific oriented functions, classes, etc, so its uses are varied. It is freely available, since it is an open source project.

ROOT comes with CINT, a C++ interpreter. This can be very useful for debugging, but unfortunately for this project, too slow, which is why to use my code it is necessary to compile it first (See section COMPILING).

The code uses a configuration file (lvis_mol_lens.conf) where all the parameters are specified. On top of that I have written a wrapper script that takes care of all the linking and facilitates the interface between BASH and ROOT. Running the code is as simple as running wrapper.C.

There are a couple of ways of running wrapper. The simplest one relies entirely on the parameters in the configuration file. In the ROOT command line do:

```
.x wrapper.C+
```

Wrapper also accepts arguments:

- arg: Acts as a "toggle" switch for telling the program to accept the following arguments. Possible values are 0 or 1.

- region: Tells the program which region of the simulation to change a parameter. Possible values are "mol" (molasses) or "lens" (laser lens).

- par_type: Parameter type to change in specified region. Possible values are: waist, detun (detuning) and w0pos (waist position).

- ite: Iteration number. This gets multiplied by step. It is useful if running in a loop.

- step: The step to vary at each iteration.

Using these arguments, one can run wrapper externally, for example, from a BASH script. This a sample script I used to run the code in batch mode (no plots) for several iterations of the detuning parameter in the molasses region:

```#!/bin/sh
for i in `seq 0 4`;
do
root.exe -b -q -l 'wrapper.C+('1',"'mol'","'detun'",'$i','0.1')'
wait
done
```

COMPILING:
This part is easy. All one has to do in a ROOT command line is:

```
.L wrapper.C+
```

It is important to keep the files read_inputs.cxx and read_inputs.h in the same folder as the rest. If changes are made to the script, it is important to check that the appropriate header files are being included in lvis_mol_lens.C.

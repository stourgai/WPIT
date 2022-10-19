# WPIT
__________________________________________________________________________________________________________________________________________________________
The Wave-Particle Interactions Toolset is an open source, Python-based set of tools for modelling the interactions between energetic charged particles and whistler-mode waves in the magnetosphere through test particle simulations. Whereas numerous modelling techniques are described in literature, there is a lack of a unified open-source toolset that incorporates the equations and parameterizations used by different wave-particle interaction models in a user-friendly environment. WPIT incorporates key routines related to wave-particle interactions in a Jupyter Notebook environment, enabling the traceability of all relevant equations in terms of their derivation and key assumptions, together with the programming environment and integrated graphics that enable users to conduct state-of-the-art wave-particle interaction simulations rapidly and efficiently. Specific capabilities of WPIT include: (a) routines for
the calculation of key parameters of the environment, such as plasma densities, frequencies and
geomagnetic field; (b) routines for the calculation of wave properties, such as the refractive index,
wave amplitude and resonance angles, based on the local environmental parameters, according
to cold plasma and fully adiabatic warm plasma theories; (c) models of the interaction of charged
particles with VLF waves under gyro-averaged equations of motion, including routines for the
characterization of the degree of nonlinearity of wave-particle interactions; and (d) a Landau
damping module, which calculates the damping of waves along the ray path due to interactions
with warm particle populations. WPIT can be used either as a stand-alone simulation tool or as
a library of routines that the user can extract and incorporate into an independent simulation.

![WPIT repository overview](wpit_overview.png)

________________________________________________________________________________________________________________________________________________________
### The WPIT Repository

The WPIT folder contains the source code builded in the form of modules and sub-modules. WPIT contains four modules, i.e. the *Environment_mod*, the *WaveProperties_mod*, the *LandauDamp_mod* and the *WPI_mod*, which includes three sub-modules, i.e. the *whistler_electron_mod*, the *EMIC_ion_mod* and the *parallel_EMIC_mod*. The *Module_descriptions* folder contains Jupyter notebooks with analytic theoretical description of the equations of each module along with example calls of each routine. The *WPI_tests* folder contains reproduction of results of the literature in Jupyter notebook formats, which act as a verification of the code and as tutorials of the use of WPIT. The *WPIT_Results* folder contains the Jupyter notebooks of the simulations presented in WPIT paper. Lastly, the *Documentation* folder includes the API documentation of the source code in .html format.

________________________________________________________________________________________________________________________________________________________

### Testing

The WPIT code has been tested in Ubuntu 18.04LST with Python 3.6.9. The version of the packages for testing are:

- matplotlib 3.2.0
- numpy 1.19.5
- scipy 1.5.4
- pandas 1.1.5
- spacepy 0.2.2
- notebook 6.4.1

________________________________________________________________________________________________________________________________________________________

### Installation

The *requirements.txt* contains all the package dependencies of WPIT. To install them, run:

```pip install -r requirements.txt```

Just clone the repository and WPIT is ready to run.

________________________________________________________________________________________________________________________________________________________

### Copyright

© 2022 Tourgaidis and Sarris.
WPIT is an open-access code distributed under the terms of the
Creative Commons Attribution License (CC BY). The use, distribution or
reproduction in other forums is permitted, provided the original
author(s) and the copyright owner(s) are credited and that the following
publication describing the code is cited, in accordance with accepted
academic practice:

Tourgaidis, S., & Sarris, T. Wave-Particle Interactions Toolset (WPIT): a Python-based toolset to model wave-particle interactions in the magnetosphere. Frontiers in Astronomy and Space Sciences, 295.


The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and non infringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.

________________________________________________________________________________________________________________________________________________________

### Funding

![Funding](hfri_funding.png)

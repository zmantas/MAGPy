'''
Written by Van Buchem in 2020

Adapted from 'MAGMA.FOR' for which the original code was written by Al Cameron in 1987 ###

History of Fortran version:
- modified by Bruce Fegley July 2002 (bfegley@wustl.edu)
- modified by Laura Schaefer August 2002, June 2003(save file) (laura_s@levee.wustl.edu)
- thermodynamic data comments added by Bruce Fegley June 2004
- modified by Laura Schaefer March 2005 - added Fe2O3 and Fe3O4,
  added commenting, and set CON(El) to zero when A(El) is zero
- modified by Yamila Miguel June 2011 (miguel@mpia-hd.mpg.de)


References describing the MAGMA code:
- B. Fegley, Jr. & A.G.W. Cameron (1987) A Vaporization Model for 
  Iron/Silicate Fractionation in the Mercury Protoplanet. 

- Earth Planet. Sci. Lett. 82, 207-222.

- L. Schaefer & B. Fegley, Jr. (2004)A Thermodynamic Model of High

- Temperature Lava Vaporization on Io. Icarus 169, 216-241.

'''

import sys
import numpy as np
from model_comp import model_comp
from simulation import simulation

def main():

    ''' Temperature for which the calculations will be done '''
    T = 1613 # Temperature of magma in Kelvin

    ''' File names '''
    input_fname = 'input/initial-composition.dat'
    output_fname = 'output/MAGMA.OUT'

    ''' Loading initial composition '''
    model = model_comp(input_fname)
    model.print()

    ''' Start simulation '''
    sim = simulation(model,T)
    sim.start()

    ''' Write to output ''' 
    sim.write_to_output(output_fname)

if __name__ == "__main__":
    sys.exit(main())


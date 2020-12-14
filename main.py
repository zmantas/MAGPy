# Standard libraries
import sys
from tqdm import tqdm

# Local modules
from melt_vapor_system import system
from melt_activity import melt_activity
from vapor_pressure import vapor_pressure
from vaporiser import vaporise


def main():

    ''' 
    Setting initial values 
    '''

    # Temperature for which the calculations will be done '''
    T = 2200 # Temperature of magma in Kelvin
    V = 0.9 

    # File names
    # input_fname = 'input/BSE-initial.dat'
    input_fname = 'input/ic_Komatiite.dat'
    output_fname = 'output/MAGMA.OUT'

    '''
    Initialising classes
    '''
    sim = system(input_fname,T)
    melt = melt_activity(sim)
    vapor = vapor_pressure(sim)
    vap = 0
    it = 0

    with tqdm(total=1) as pbar:
        pbar.set_description(f'Vaporization percentage (stops at {V*100}):')

        while vap < V and it <= 1e5:

            # Calculating activities and partial pressures
            melt.melt_activity_calculation(sim)
            vapor.vapor_pressure_calculation(sim)

            # Reapeating calculations by adding F2O3
            melt.melt_activity_calculation(sim,addF2O3=True)
            vapor.vapor_pressure_calculation(sim)

            # TODO: Output first equilibrium before removal of vapor

            # Remove vapor
            vap = vaporise(sim,vapor)

            # Update counters
            it += 1
            pbar.update(vap - pbar.n)

    pbar.close()

if __name__ == "__main__":
    sys.exit(main())
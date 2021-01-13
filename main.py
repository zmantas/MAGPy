# Standard libraries
import sys
from tqdm import tqdm

# Local modules
from library.melt_vapor_system import system
from library.melt_activity import melt_activity
from library.vapor_pressure import vapor_pressure
from library.vaporiser import vaporise
import library.print_functions as print_functions
import magpy_cfg

def main():

    # Setting initial values 
    T = magpy_cfg.magmaT # Temperature of magma in Kelvin
    V = magpy_cfg.vaporFrac  # Set desired vaporisation fraction (must be < 1)

    # File names
    if magpy_cfg.comp == 'BSE':
        input_fname = 'input/BSE.dat'
    elif magpy_cfg.comp == 'Komatiite':
        input_fname = 'input/Komatiite.dat'
        
        
    output_fname = 'output/magpy.out'
    outputEle_fname = 'output/magpyVapor.out'

    # Initialising classes and trackers
    sim = system(input_fname,T)
    melt = melt_activity(sim)
    vapor = vapor_pressure(sim)
    vap = 0
    it = 0

    # Printing initial parameters
    print_functions.print_init(sim,output_fname)

    # tqdm progres bar
    with tqdm(total=1) as pbar:
        pbar.set_description(f'Vaporization percentage (stops at {int(V*100)}%)')

        while vap < V and it <= 1e5 or it == 0:

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

    print_functions.print_results(sim,melt,vap,output_fname)
    print_functions.print_resultsEle(sim,melt,vap,outputEle_fname)

if __name__ == "__main__":
    sys.exit(main())
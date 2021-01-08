# Standard libraries
import sys
import numpy as np
from tqdm import tqdm

# Local modules
from library.melt_vapor_system import system
from library.melt_activity import melt_activity
from library.vapor_pressure import vapor_pressure
from library.vaporiser import vaporise
import library.print_functions as print_functions

def calc(T,V,input_fname,output_fname,file,first=False):

    # Initialising classes and trackers
    sim = system(input_fname,T)
    melt = melt_activity(sim)
    vapor = vapor_pressure(sim)
    vap = 0
    it = 0

    # Printing initial parameters
    print_functions.print_init(sim,output_fname)

    # tqdm progres bar
    # with tqdm(total=1) as pbar:
    #     pbar.set_description(f'Vaporization percentage (stops at {int(V*100)}%)')

    while vap < V and it <= 1e5 or it == 0:

        # Calculating activities and partial pressures
        melt.melt_activity_calculation(sim)
        vapor.vapor_pressure_calculation(sim)

        # Repeating calculations by adding F2O3
        melt.melt_activity_calculation(sim,addF2O3=True)
        vapor.vapor_pressure_calculation(sim)

        # TODO: Output first equilibrium before removal of vapor

        # Remove vapor
        vap = vaporise(sim,vapor)

        # Update counters
        it += 1
        # pbar.update(vap - pbar.n)

    # pbar.close()

    if first:
        for ox in sim.gasMoleFrac:
            file.write(f'{ox},')
        file.write('\n')

    for ox in sim.gasMoleFrac:
        file.write(f'{sim.gasMoleFrac[ox]},')
    file.write('\n')


    #print_functions.print_results(sim,melt,vap,output_fname)

def main():

    # Setting initial values 
    T = [1500,3000] # Temperature range of magma in Kelvin
    V = 0  # Set desired vaporisation fraction (0 <= V < 1)

    print(f'Vaporization percentage set to: {V}')
    print(f'Calculating for temperature range: {T[0]}-{T[1]} K')

    # File names
    input_fname = 'input/BSE-initial.dat'
    output_fname = 'output/MAGMA.OUT'

    # Make output file for temp var
    file = open('output/temp_var.csv', 'w')

    # Run calculations for entire temperature range
    for t in tqdm(np.arange(T[0],T[1]+1)):
        if t == T[0]:
            calc(t,V,input_fname,output_fname,file,first=True)            
        calc(t,V,input_fname,output_fname,file)

    file.close()
    

if __name__ == "__main__":
    sys.exit(main())
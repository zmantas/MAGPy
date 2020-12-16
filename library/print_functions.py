

def print_init(sim,output_fname):

    ''' Opening the file '''
    print(f'Writing output to: {output_fname}')
    file = open(output_fname, "w")

    file.write(f'INITIAL PARAMETERS AND COMPOSITION\n\n')

    ''' Magma composition ''' 
    file.write(f'Magma composition:\n \n')
    file.write(f'Oxide        WT%          Mole%\n')
    for ox in sim.comp_init:
        file.write(f'{ox:<5}    {sim.comp_init[ox]:<9.4e}    {sim.molOx_perc[ox]:.6e}\n') 
    file.write(f'Total    {sim.totWt:<9.4e}    {sim.totPerc:.6e}\n')

    ''' Abundances '''
    file.write('\nAtomic abundances on cosmochemical scale\n\n')
    for el in sim.abEl:
        file.write(f'{el:<4}=  {sim.abEl[el]:.6e}\n')
    ''' Oxide Mole fraction in silicate '''
    file.write(f'\nOxide Mole Fraction (F) in Silicate\n \n')
    for name in sim.fAbOx:
        file.write(f'{sim._metal2oxide[name]:<6}= {sim.fAbOx[name]:.6e}\n') 

    ''' Relative atomic abundances of metals '''
    file.write(f'\nRelative atomic abundances of metals\n \n')
    for name in sim.fAbAtom:
        file.write(f'{name:<3}= {sim.fAbAtom[name]:.6e}\n')

    file.close()

def print_results(sim,melt,vap,output_fname):

    file = open(output_fname, "a")

    file.write(f'\nFINAL COMPOSITION\n')

    ''' The activity data '''
    file.write(f'\nActivity coefficients (G) of oxides in the melt\n \n')
    for ox in sim.gamma:
        file.write(f'{ox:<4}= {sim.gamma[ox]:.6e}\n')

    file.write(f'\nActivities (A) of Species in the melt\n \n')
    for ox in sim.act_ox:
        file.write(f'{ox:<12} {sim.act_ox[ox]:.6e}\n')
    for spec in melt.act_pseudo:
        file.write(f'{melt.name_pseudo[spec]:<12} {melt.act_pseudo[spec]:.6e}\n') 

    ''' Gas pressure partial vapor ''' 
    file.write('\nGas partial pressures (P) in vapor \n\n')
    file.write(f'{"O":<8} {sim.presGas["O"]:.6e}\n')
    file.write(f'{"O2":<8} {sim.presGas["O2"]:.6e}\n')
    for name_el in sim.abEl:  
        for gas in sim.presGas:
            if name_el in gas:
                file.write(f'{gas:<8} {sim.presGas[gas]:.6e}\n')
    file.write(f'{"e-":<8} {sim.presGas["EnE"]:.6e}\n')
    file.write(f'\n{"Total":<8} {sim.totPres:.6e}\n')

    ''' Gas mole fractions (X) in vapor ''' 
    file.write('\nGas mole fractions (X) in vapor \n\n')
    file.write(f'{"O":<8} {sim.gasMoleFrac["O"]:.6e}\n')
    file.write(f'{"O2":<8} {sim.gasMoleFrac["O2"]:.6e}\n')
    for name_el in sim.abEl:  
        for gas in sim.gasMoleFrac:
            if name_el in gas:
                file.write(f'{gas:<8} {sim.gasMoleFrac[gas]:.6e}\n')
    file.write(f'{"e-":<8} {sim.gasMoleFrac["EnE"]:.6e}\n')

    ''' Vaporized fractions '''
    file.write('\nFraction of magma that is vaporized\n\n')
    file.write(f'Vap. fraction  = {vap:.6e}\n')
    file.write(f'Weight percent = {sim.massFrac:.6e}\n')

    file.close()
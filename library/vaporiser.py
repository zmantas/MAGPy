
def vaporise(sim,vapor):

    # Calculate gas pressure for each element
    tot_gasPres = {}
    for name_el in sim.abEl:
        tot_gasPres[name_el] = 0
        # Selects al gasses that contain the relevant element
        for name_gas in sim.presGas:
            if name_el in name_gas:
                tot_gasPres[name_el] += sim.presGas[name_gas]

    # Added seperatly due to absence from original abundance list
    tot_gasPres['O'] = sim.presGas['O'] + sim.presGas['O2']  
    tot_gasPres['EnE'] = sim.presGas['EnE'].copy() 
    sim.totPres = sum(tot_gasPres.values()) 

    # Gas mole fraction (P/Ptot)
    sim.gasMoleFrac = {gas : sim.presGas[gas]/sim.totPres for gas in sim.presGas}
    
    # Calculating total mole fractions of each element in the gas
    n_tot = sum(vapor.n_el.values())
    n_frac = {element : vapor.n_el[element]/n_tot for element in vapor.n_el}

    # Calculating volatilities using mole fraction and atomic abundance
    volatilities = {element: n_frac[element]/sim.fAbAtom[element] \
                    if sim.fAbAtom[element] > 1e-20 else 0 for element in sim.fAbAtom}

    # Calculating vaporisation fraction using most volatile element
    vapoFrac = 0.05/max(volatilities.values())

    # Calculating volatilities using mole fraction and elemental abundance
    volatilities1 = {element: n_frac[element]/sim.abEl[element] \
                     if sim.fAbAtom[element] > 1e-20 else 0 for element in sim.abEl}
    
    # Calculating vaporisation fraction using most volatile element
    vapoFrac1 = 0.05/max(volatilities1.values())        

    # Calculate new abundances for each element 
    sim.fAbAtom = {element: sim.fAbAtom[element] - vapoFrac * n_frac[element] \
                        for element in sim.fAbAtom}

    sim.abEl = {element: sim.abEl[element] - vapoFrac1 * n_frac[element]\
                     for element in sim.abEl}

    # Ensuring that no abundance values are negative
    for element in sim.abEl:
        if element == 'Fe':
            if sim.abEl[element] <= 2e-20:
                sim.abEl[element] = 0
                sim.fAbAtom[element] = 0    
        else:
            if sim.abEl[element] <= 0:
                sim.abEl[element] = 0
                sim.fAbAtom[element] = 0

    # Fraction of vaporised magma (new_totAbundance/old_totAbundance)
    vap = 1 - sum(sim.abEl.values())/sim.abETot 

    # Calculate weight percent vaporized 
    massVapo = 0
    for element in sim.abEl:
        for oxide in sim._mwOxides:
            if element in oxide:
                if element in ['Al','Na','K']:
                    massVapo += 0.5 * sim.abEl[element] * sim._mwOxides[oxide]
                elif oxide != 'Fe2O3':
                    massVapo += sim.abEl[element] * sim._mwOxides[oxide]
    sim.massFrac = (sim.mass-massVapo)/sim.mass

    # Renormalize the abundances 
    fAbAtom_tot = sum(sim.fAbAtom.values()) # Mole fraction of elemens
    fAbOx_tot = 0  # Mole fraction of oxides
    for el in sim.fAbAtom:
        if el in ['Al','Na','K']:
            fAbOx_tot += 0.5 * sim.fAbAtom[el]
        else:
            fAbOx_tot += sim.fAbAtom[el] 

    for el in sim.fAbOx: # Relative abundnaces of metals by moleculre
        if el in ['Al','Na','K']:
             sim.fAbOx[el] = 0.5 * sim.fAbAtom[el]/fAbOx_tot
        else:
             sim.fAbOx[el] = sim.fAbAtom[el]/fAbOx_tot 

    # Relative abundances of metals by atom
    sim.fAbAtom = {el : sim.fAbAtom[el]/fAbAtom_tot for el in sim.fAbAtom}

    return vap

import numpy as np
from data.thermodynamic_data import thermodynamic_data

class simulation():
    ''' 
    Initialises the vaporisation simulation.
    '''
    def __init__(self,model_composition,T):

        ''' File names ''' 
        mwOxides_fname = 'data/weights_oxides.csv'
        wMetals_fname = 'data/weights_metals.csv'

        ''' Setting temperature '''
        self.T = T   # temperature of magma (Kelvin)
        print(f'\nCalculating for a magma temperature of {self.T} K')
        
        ''' Importing metal and oxide name lists from '''
        self._metalNames = model_composition._metalNames
        self._oxideNames = model_composition._oxideNames

        ''' Relating oxides to metals '''
        self._metal2oxide = {}
        for i,metal in enumerate(self._metalNames):
            self._metal2oxide[metal] = self._oxideNames[i]

        '''Constants'''
        self._avog = 6.023e23 # Avogadro's number

        '''Setting simulation parameters TODO: Make variable'''
        self.iRep = 0  
        self.iStep = 0  # counts the number of vaporization steps
        self.iPrn = 0   # used to determine which steps are printed to the output file
        self.iIt = 0    # counts the number of iterations needed to solve for the activities
        self.iFirst = 0
        self.perVap = 0 # number of iterations of vaporization process

        '''Importing molecular oxide and metal weights'''
        self._mwOxides = dict(np.genfromtxt(mwOxides_fname,delimiter= ',',names=True,skip_header=1,\
                                 dtype=None, encoding = 'UTF-8'))
        self._wMetals = dict(np.genfromtxt(wMetals_fname,delimiter= ',',names=True,skip_header=1,\
                                 dtype=None, encoding = 'UTF-8',autostrip=True))

        # Evaluating the multiplications and dividing by _avog (mol wt. in g/mole)
        self._wMetals.update((x, eval(y)/self._avog) for x, y in self._wMetals.items()) 

        ''' Calculating number of moles of the oxides calculated from the input composition '''
        
        # Total weight
        self.wt_oxides = model_composition.wt_oxides
        self.totWt = sum(self.wt_oxides.values())
        
        # Number of moles per oxide
        self.numMolOx = {}
        for ox in self._mwOxides:
            if self.wt_oxides[ox] != 0: # Avoiding div by 0 error
                self.numMolOx[ox] = self.wt_oxides[ox] / self._mwOxides[ox]
            else: 
                self.numMolOx[ox] = 0
                
        # Total amount of moles
        self.totMol = sum(self.numMolOx.values())
        
        # Mole percentages of the oxides
        self.mPerc = {}
        for ox in self.numMolOx:
            if self.wt_oxides[ox] != 0: # Avoiding div by 0 error
                self.mPerc[ox] = 100 * self.numMolOx[ox] / self.totMol
            else:
                self.mPerc[ox] = 0

        # Total percentage
        self.totPerc = sum(self.mPerc.values())

        ''' 
        If there is Si in the composition, abundances are normalised to Si = 1e6 
        AE(element) = (# of metal atoms in oxide) * MO(moles of oxide)
        * 1D6 / moles of SiO2
        '''
        self.abEl = {} # abundances of the elements

        if self.numMolOx['SiO2'] != 0:            
            self.abEl['Si'] = 1e6
            self.abEl['Mg'] = self.numMolOx['MgO']  * 1e6 / self.numMolOx['SiO2']
            self.abEl['Al'] = self.numMolOx['Al2O3'] * 2 * 1e6 / self.numMolOx['SiO2']
            self.abEl['Ti'] = self.numMolOx['TiO2'] * 1e6 / self.numMolOx['SiO2']
            self.abEl['Fe'] = (self.numMolOx['FeO'] + 2 * self.numMolOx['Fe2O3']) * 1e6 / self.numMolOx['SiO2']
            self.abEl['Ca'] = self.numMolOx['CaO']  * 1e6 / self.numMolOx['SiO2']
            self.abEl['Na'] = self.numMolOx['Na2O'] * 2 * 1e6 / self.numMolOx['SiO2']
            self.abEl['K']  = self.numMolOx['K2O']  * 2 * 1e6 / self.numMolOx['SiO2']

        else:
            self.abEl['Si'] = self.numMolOx['SiO2'] * self._avog
            self.abEl['Mg'] = self.numMolOx['MgO']  * self._avog
            self.abEl['Al'] = self.numMolOx['Al2O3'] * 2 * self._avog
            self.abEl['Ti'] = self.numMolOx['TiO2'] * self._avog
            self.abEl['Fe'] = (self.numMolOx['FeO'] + 2 * self.numMolOx['Fe2O3']) * self._avog
            self.abEl['Ca'] = self.numMolOx['CaO']  * self._avog
            self.abEl['Na'] = self.numMolOx['Na2O'] * 2 * self._avog
            self.abEl['K']  = self.numMolOx['K2O']  * 2 * self._avog  

        ''' Initial weight of melt so that weight% vaporized can be calculated later '''
        self.mass = 0
        for element in self.abEl:
            for oxide in self._mwOxides:
                if element in oxide:
                    if element in ['Al','Na','K']:
                        self.mass += 0.5 * self.abEl[element] * self._mwOxides[oxide]
                    elif oxide != 'Fe2O3':
                        self.mass += self.abEl[element] * self._mwOxides[oxide]

        ''' Renomarlizing the abundances ''' 
        # Total atomic abundance of all the elements (except 0)
        self.abETot = sum(self.abEl.values())
        # Molecular abundance of all the oxides
        self.abETot1 = self.abEl['Si'] + \
                  self.abEl['Mg'] + \
                  self.abEl['Fe'] + \
                  self.abEl['Ca'] + \
                  0.5 * self.abEl['Al'] + \
                  self.abEl['Ti'] + \
                  0.5 * self.abEl['Na'] + \
                  0.5 * self.abEl['K'] # May have to check if halving these values is ok

        ''' 
        Relative abundances of the metals (by molecule and atom) in the mantle 
        (same as oxide mole fractions in the magma)
        ''' 
        self.fAbMolecule = {} # Relative abundance of the metals per molecule
        self.fAbAtom = {}     # Relative abundance of the metals per atom

        for ox in self.abEl:
            if ox in ['Al','Na','K']: # Halves values for Al, Na and K
                self.fAbMolecule [ox] = 0.5 * self.abEl[ox] / self.abETot1
            else:   
                self.fAbMolecule [ox] = self.abEl[ox] / self.abETot1
            self.fAbAtom[ox] = self.abEl[ox] / self.abETot

        ''' Setting initial values of key pressures and adjustment factors for gases '''
        self.presGas = {}  # gas pressures (initialised below)
        self.presLiq = {} # liquid pressures (not defined until start of calc)
        self.adjFact = {}  # adjustment factors
        self.gas_names = ['SiO','O2','MgO','Fe','Ca','Al','Ti','Na','K']
        for gas in self.gas_names:
            self.presGas[gas] = 1
            self.adjFact[gas] = 1
        
        ''' Setting inital values of gamma (activity coefficient)'''
        self.gamma = {}  # activity coefficient - gamma
        self.liquid_names = ['Si','Mg','Fe','Fe3','Ca','Al','Ti','Na','K'] # TODO: Check name
        for liquid in self.liquid_names: 
            self.gamma[liquid] = 1 # TODO: May have to do same for second parameter (l318 in og code)

        ''' Computing activities in the melt (based on temp) '''
        self.td = thermodynamic_data(self.T)
        self.actMelt = self.td.activities_melt()

    # end __init__()

    def activity_gas_calculation(self,addF2O3=True): 
        # print(addF2O3)
        self.iit = 0 # Counter needed to solve for activities
        while self.iit < 1e8: # Setting to 1e5 purely so that no unending loop is created, check best size of value 
            print('IIT',self.iit, self.iRep, self.iStep)
            ''' 
            Activities of oxides in the melt:
            - formula: activity = molecular abundance * activity coefficient
            '''
            # if self.iStep == 806:
            #     print('Before actox calc')
            #     print(self.fAbMolecule['Na'], self.gamma['Na'])

            self.actOx = {}
            for i,metal in enumerate(self._metalNames):
                if metal != 'Fe3':
                    self.actOx[self._oxideNames[i]] = self.fAbMolecule[metal] * self.gamma[metal]
                else:
                # Activity of Fe2O3 is estimated using gas chemistry, then all acitivities are recomputed
                    if addF2O3:
                        self.actOx['Fe2O3'] = self.presLiq['Fe2O3'] * self.gamma['Fe3']
                        # print('Fe2O3',self.actOx['Fe2O3'])
                        # print(self.presLiq['Fe2O3'],self.gamma['Fe3'])
                    else:
                        self.actOx['Fe2O3'] = 0
            # print(f'SiO2 {self.actOx["SiO2"]}, MgO {self.actOx["MgO"]}')
            # print(f'Fe2O3 {self.actOx["Fe2O3"]}')

            ''' Computes activities of the complex melts '''        
            self.actMeltComp = self.td.activities_meltComplex(self.actOx,self.iStep)

            ''' Recomputing gamma from the activities computed above '''
            self.gamma_new = self.td.recompute_gamma(self.actOx,self.gamma,addF2O3,self.fAbMolecule,self.iStep)
            
            ''' Compute ratio of newly computed activity and previous activity '''
            self.gamRat = {}
            for element in self.gamma:
                if self.gamma_new[element] != 0: # TODO: Check if for Fe3 you need to check that Fe is zero as well
                    self.gamRat[element] = self.gamma_new[element]/self.gamma[element]
                else:
                    self.gamRat[element] = 1

            ''' 
            If self.gamRat[elements] ~1, the code has arrived at a solution for all the 
            activities, and moves on to the gas chemistry. 
            If this is not the case, then the activity coefficients are adjusted and 
            the activities are recomputed until a solution is foudn. 
            ''' 
            # print(list(self.gamRat.values()))
            if np.sum(np.isnan(list(self.gamRat.values())))>0:
                print('gamrat',self.gamRat)
                print('gamnew',self.gamma_new)
                print('gamma',self.gamma)
                raise RuntimeError('WTF BRO THERE IS NAN BREAD')

            if all(rat < 1e-5 for rat in np.abs(np.log10(list(self.gamRat.values())))):
                #print(f'Solution for the gas activities found after {self.iit+1} iteration(s).')
                break

            if self.iit > 500:
                for element in self.gamma_new:
                    self.gamma_new[element] = (self.gamma_new[element] * self.gamma[element]**4)**(1/5)
            elif self.iit > 30:
                for element in self.gamma_new:
                    self.gamma_new[element] = (self.gamma_new[element] * self.gamma[element]**2)**(1/3)
            else:
                for element in self.gamma_new:
                    self.gamma_new[element] = (self.gamma_new[element] * self.gamma[element])**(1/2)

            self.gamma = self.gamma_new.copy() # Update gamma values
            self.iit += 1 # updating counter
            
            # print(f'gamma Si {self.gamma["Si"]:.6e} Mg {self.gamma["Mg"]:.6e} Fe Si {self.gamma["Fe"]:.6e}')
            # print(f'iit {self.iit}')
            if self.iit >= 1e8: 
                raise RuntimeError('Max recursion limit reached while calculating activities.')


        # end while loop

        '''
        .....................................................................
        GAS CHEMISTRY CALCULATIONS
        Calculate gas chemistry in equilibrium with the calculated activities
        from above.
        .....................................................................

        The adjustment factors are recomputed and the gas chemistry is repeated until 
        a solution where the adjFact ~1 or 0 is converged upon.
        '''
        self.iit = 0 # counter
        self.dif_range = 2.30359e-6 # Max amount that the adjFact can difer from 1
        while not all( 1-self.dif_range < fact < 1+self.dif_range or fact == 0 \
                   for fact in list(self.adjFact.values())) or self.iit == 0:
            '''
            ADJUST THE ABUNDANCES OF THE MAJOR GASES OF EACH ELEMENT
            these abundances are used to calculate all other gas chemistry
            '''
            for gas in self.gas_names:
                self.presGas[gas] = self.presGas[gas]* self.adjFact[gas]

            ''' Compute the partial pressures of the vapor species ''' 
            self.td.ion_chemistry(self.presGas,self.presLiq)

            ''' Calculate the number densities of each species and for each element'''
            self.td.number_density(self.presGas)

            ''' Recompute the adjustment factors for the key pressures ''' 
            self.td.recompute_adjFact(self.presGas,self.presLiq,self.gamma,self.adjFact,self.fAbMolecule,self.actOx)
            
            ''' While loop break '''
            self.iit += 1 # updating counter
            if self.iit >= 1e8: 
                raise RuntimeError('Max recursion limit reached while calculating adjustment factors.')

        #print(f'Solution for the adjustment factors found after {self.iit} iteration(s).')

        # end while loop

    # end activity_gas_calculation

    def vaporisation(self):

        ''' Calculating total gas pressure '''
        self.tot_gasPres = {}
        for name_el in self.abEl:
            self.tot_gasPres[name_el] = 0
            for name_gas in self.presGas:
                if name_el in name_gas:
                    self.tot_gasPres[name_el] += self.presGas[name_gas]

        self.tot_gasPres['O'] = self.presGas['O'] + self.presGas['O2'] # Adding the oxygen pressure 
        self.tot_gasPres['EnE'] = self.presGas['EnE'].copy() 
        self.totPres = sum(self.tot_gasPres.values()) # total pressure

        ''' Calculating the gas mole fractions (P/PTOT) '''
        self.gasMoleFrac = {}
        for gas in self.presGas:
            self.gasMoleFrac[gas] = self.presGas[gas]/self.totPres

        ''' Calculate the total mole fraction of an element in the gas '''
        self.sum_totRho = sum(self.td.totRho.values())
        self.totMolFracEl = {}
        for element in self.td.totRho:
            self.totMolFracEl[element] = self.td.totRho[element]/self.sum_totRho

        # if self.addF2O3 == 1:
            #TODO: print the calculated aquilibrium abundances before removing mass for the
            #      first vaporisation step
            # pass

        '''
        Compute the step size.  The most volatile element in the melt will be 
        reduced by 5%.  The computation is done for both the mole fractions of
        the elements and the atomic abundances because once PLAN(element) (self.abEl) becomes
        smaller than ~1D-35, it gets put to zero due to the limitations of FORTRAN.
        The mole fractions are renormalized below and are used to compute the 
        oxide mole fractions.  The PLAN(element) values are 
        used to compute the fraction vaporized.
        TODO: CHECK IF THIS STILL HOLDS FOR PYTHON
        '''

        # Calculating volatilities using mole fraction and atomic abundance
        volatilities = {element: self.totMolFracEl[element]/self.fAbAtom[element] if self.fAbAtom[element] > 1e-20 else 0 for element in self.fAbAtom}
        # Calculating vaporisation fraction using most volatile element
        # print(volatilities.values())
        # print(max(volatilities.values()))
        self.vapoFrac = 0.05/max(volatilities.values())


        # Calculating volatilities using mole fraction and elemental abundance
        volatilities1 = {element: self.totMolFracEl[element]/self.abEl[element] if self.fAbAtom[element] > 1e-20 else 0 for element in self.abEl}
        # Calculating vaporisation fraction using most volatile element
        # print(volatilities.values())
        # print(max(volatilities.values()))
        self.vapoFrac1 = 0.05/max(volatilities1.values())        

        '''
        Compute the new abundance of each element to be used in the next
        vaporization step. 
        CON(El) = mole fraction in melt
        Plan(El) = total atomic abundance
        if plan(el) = 0, then set con(el) = 0 also.
        '''
        self.fAbAtom = {element: self.fAbAtom[element] - self.vapoFrac*self.totMolFracEl[element] \
                        for element in self.fAbAtom}

        self.abEl = {element: self.abEl[element] - self.vapoFrac1 * self.totMolFracEl[element]\
                     for element in self.abEl}

        for element in self.abEl:
            if element == 'Fe':
                if self.abEl[element] <= 2e-20:
                    self.abEl[element] = 0
                    self.fAbAtom[element] = 0    
            else:
                if self.abEl[element] <= 0:
                    # print(element, ' less than 0')
                    self.abEl[element] = 0
                    self.fAbAtom[element] = 0
                    # print('check it',self.abEl[element],self.fAbAtom[element])

        self.abETot_new = sum(self.abEl.values())
        self.abRatio = self.abETot_new/self.abETot
        self.vap = 1 - self.abRatio # fraction of magma that is vaporized
        
        ''' Calculate weight percent vaporized '''
        self.massVapo = 0
        for element in self.abEl:
            for oxide in self._mwOxides:
                if element in oxide:
                    if element in ['Al','Na','K']:
                        self.massVapo += 0.5 * self.abEl[element] * self._mwOxides[oxide]
                    elif oxide != 'Fe2O3':
                        self.massVapo += self.abEl[element] * self._mwOxides[oxide]
        self.massFrac = (self.mass-self.massVapo)/self.mass

        ''' Renormalize the abundances '''
        self.conTot_el = sum(self.fAbAtom.values()) # Mole fraction of elemens
        self.conTot_ox = 0  # Mole fraction of oxides
        for el in self.fAbAtom:
            if el in ['Al','Na','K']:
                 self.conTot_ox += 0.5 * self.fAbAtom[el]
            else:
                self.conTot_ox += self.fAbAtom[el] 

        for el in self.fAbMolecule: # Relative abundnaces of metals by moleculre
            if el in ['Al','Na','K']:
                 self.fAbMolecule[el] = 0.5 * self.fAbAtom[el]/self.conTot_ox
            else:
                 self.fAbMolecule[el] = self.fAbAtom[el]/self.conTot_ox 
        # print(self.fAbMolecule['Na'], self.abEl['Na'])
        # print(self.fAbAtom['Na'],self.conTot_ox)
        # Relative abundances of metals by atom
        self.fAbAtom = {el : self.fAbAtom[el]/self.conTot_el for el in self.fAbAtom}

    # end vaporisation

    def start(self):

        self.vap = 0

        '''
        If this is the first run through activity/gas calculations for this step
        then go back to activity calculations and add in the iron oxides
        Fe2O3 and Fe3O4. Repeat activity and gas calculations until answers converge. 
        '''
        
        file = open('magpy_trackers.csv', "w")

        while self.vap < 1 and self.iRep <= 2647:
            # print(self.iRep)
            self.activity_gas_calculation(addF2O3=False)
            self.activity_gas_calculation()
            # print(self.gamma['Si']) 
            self.vaporisation()

            # print(f'{self.vap:.6e}')

            file.write(f'{self.iRep} {self.vap:.6e}\n')

            ''' Update counters '''
            self.iPrn += 1 # Used to track how often to print
            self.iStep += 1 # Used to track nuber of vaporization steps
            self.iRep += 1 

        # CHECK IF self.vap => 1
        # IF NOT go back to calculating activities

    # end start()

    def set_temp(self,T):
        self.T = T
        print(f'Temperature set to {self.T} (K)')

    # end set_temp()

    def set_perVap(self,perVap):
        self.perVap = perVap
        print(f'Number of vaporization permutations set to {self.perVap}')

    # end set_perVap()

    def print_init(self,output_fname):

        ''' Opening the file '''
        print(f'Writing output to: {output_fname}')
        file = open(output_fname, "w")

        ''' Printing output to file ''' 
        file.write(f'INITIAL PARAMETERS AND COMPOSITION\n')

        ''' Magma composition ''' 
        file.write(f'Magma composition:\n \n')
        file.write(f'Oxide        WT%          Mole%\n')
        for ox in self.wt_oxides:
            file.write(f'{ox:<5}    {self.wt_oxides[ox]:<9.4e}    {self.mPerc[ox]:.6e}\n') 
        file.write(f'Total    {self.totWt:<9.4e}    {self.totPerc:.6e}\n')

        ''' Abundances '''
        file.write('\nAtomic abundances on cosmochemical scale\n\n')
        for el in self.abEl:
            file.write(f'{el:<4}=  {self.abEl[el]:.6e}\n')
        ''' Oxide Mole fraction in silicate '''
        file.write(f'\nOxide Mole Fraction (F) in Silicate\n \n')
        for name in self.fAbMolecule:
            file.write(f'{self._metal2oxide[name]:<6}= {self.fAbMolecule[name]:.6e}\n') 

        ''' Relative atomic abundances of metals '''
        file.write(f'\nRelative atomic abundances of metals\n \n')
        for name in self.fAbAtom:
            file.write(f'{name:<3}= {self.fAbAtom[name]:.6e}\n')

        file.close()

    def print_results(self,output_fname):

        file = open(output_fname, "a")

        ''' The activity data '''
        file.write(f'\nActivity coefficients (G) of oxides in the melt\n \n')
        for ox in self.gamma:
            file.write(f'{ox:<4}= {self.gamma[ox]:.6e}\n')

        file.write(f'\nActivities (A) of Species in the melt\n \n')
        for ox in self.actOx:
            file.write(f'{ox:<12} {self.actOx[ox]:.6e}\n')
        for comp in self.actMeltComp:
            file.write(f'{self.td.nameMeltComp[comp]:<12} {self.actMeltComp[comp]:.6e}\n') 

        ''' Gas pressure partial vapor ''' 
        file.write('\nGas partial pressures (P) in vapor \n\n')
        file.write(f'{"O":<8} {self.presGas["O"]:.6e}\n')
        file.write(f'{"O2":<8} {self.presGas["O2"]:.6e}\n')
        for name_el in self.abEl:  
            for gas in self.presGas:
                if name_el in gas:
                    file.write(f'{gas:<8} {self.presGas[gas]:.6e}\n')
        file.write(f'{"e-":<8} {self.presGas["EnE"]:.6e}\n')
        file.write(f'\n{"Total":<8} {self.totPres:.6e}\n')

        ''' Gas mole fractions (X) in vapor ''' 
        file.write('\nGas mole fractions (X) in vapor \n\n')
        file.write(f'{"O":<8} {self.gasMoleFrac["O"]:.6e}\n')
        file.write(f'{"O2":<8} {self.gasMoleFrac["O2"]:.6e}\n')
        for name_el in self.abEl:  
            for gas in self.gasMoleFrac:
                if name_el in gas:
                    file.write(f'{gas:<8} {self.gasMoleFrac[gas]:.6e}\n')
        file.write(f'{"e-":<8} {self.gasMoleFrac["EnE"]:.6e}\n')

        ''' Vaporized fractions '''
        file.write('\nFraction of magma that is vaporized\n\n')
        file.write(f'Vap. fraction  = {self.vap:.6e}\n')
        file.write(f'Weight percent = {self.massFrac:.6e}\n')

        file.close()

    # end write_to_output()

# end simulations()

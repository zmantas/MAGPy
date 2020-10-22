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
        self.iRep = 2   
        self.iStep = 0  # counts the number of vaporization steps
        self.iPrn = 0   # used to determine which steps are printed to the output file
        self.iIt = 0    # counts the number of iterations needed to solve for the activities
        self.iFirst = 0
        self.perVap = 0 # number of iterations of vaporization process
        self.count = 0 # counting factor

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

    def activity_gas_calculation(self,count): 

        self.iit = 0 # Counter needed to solve for activities
        while self.iit < 1e5: # Setting to 1e5 purely so that no unending loop is created, check best size of value 
            ''' 
            Activities of oxides in the melt:
            - formula: activity = molecular abundance * activity coefficient
            '''
            self.actOx = {}
            for i,metal in enumerate(self._metalNames):
                if metal != 'Fe3':
                    self.actOx[self._oxideNames[i]] = self.fAbMolecule[metal] * self.gamma[metal]
                else:
                # Activity of Fe2O3 is estimated using gas chemistry, then all acitivities are recomputed
                    if count == 1:
                        self.actOx['Fe2O3'] = 0
                    else:
                        self.actOx['Fe2O3'] = self.presLiq['Fe2O3'] * self.gamma['Fe3']

            ''' Computes activities of the complex melts '''        
            self.actMeltComp = self.td.activities_meltComplex(self.actOx)

            ''' Recomputing gamma grom the activities computed above '''
            self.gamma_new = self.td.recompute_gamma(self.actOx,self.gamma,self.count,self.fAbMolecule)
            
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
            if all(rat < 1e-5 for rat in np.abs(np.log10(list(self.gamRat.values())))):
                print(f'Solution for the gas activities found after {self.iit+1} iteration(s).')
                break

            if self.iit > 50:
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
            if self.iit >= 1e5: 
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
            if self.iit >= 1e5: 
                raise RuntimeError('Max recursion limit reached while calculating adjustment factors.')

        print(f'Solution for the adjustment factors found after {self.iit} iteration(s).')

        # end while loop

    # end activity_gas_calculation

    def start(self):

        '''
        If this is the first run through activity/gas calculations for this step
        then go back to activity calculations and add in the iron oxides
        Fe2O3 and Fe3O4. Repeat activity and gas calculations until answers converge. 
        '''
        self.count += 1 # Update counter 
        self.activity_gas_calculation(self.count) 
        self.count += 1
        print('\nAdding iron oxides and repeating calculations.')
        self.activity_gas_calculation(self.count)
        

    # end start()

    def set_temp(self,T):
        self.T = T
        print(f'Temperature set to {self.T} (K)')

    # end set_temp()

    def set_perVap(self,perVap):
        self.perVap = perVap
        print(f'Number of vaporization permutations set to {self.perVap}')

    # end set_perVap()

    def write_to_output(self,output_fname):

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

        ''' Oxide Mole fraction in silicate '''
        file.write(f'\nOxide Mole Fraction (F) in Silicate\n \n')
        for name in self.fAbMolecule:
            file.write(f'{self._metal2oxide[name]:<6}= {self.fAbMolecule[name]:.6e}\n') 

        ''' Relative atomic abundances of metals '''
        file.write(f'\nRelative atomic abundances of metals\n \n')
        for name in self.fAbAtom:
            file.write(f'{name:<3}= {self.fAbAtom[name]:.6e}\n')

        ''' The activity data '''
        file.write(f'\nActivity coefficients (G) of oxides in the melt\n \n')
        for ox in self.gamma:
            file.write(f'{ox:<4}= {self.gamma[ox]:.6e}\n')

        file.write(f'\nActivities (A) of Species in the melt\n \n')
        for ox in self.actOx:
            file.write(f'{ox:<12} {self.actOx[ox]:.6e}\n')
        for comp in self.actMeltComp:
            file.write(f'{self.td.nameMeltComp[comp]:<12} {self.actMeltComp[comp]:.6e}\n') 

        ''' Closing the file ''' 
        file.close()

    # end write_to_output()

# end simulations()

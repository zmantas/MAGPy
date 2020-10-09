import numpy as np

class simulation():
    ''' 
    Initialises the vaporisation simulation.
    '''
    def __init__(self):

        '''Constants'''
        self._pConv = 1.01325e6/1.38046e-16 # _pConv converts the pressures into number densities
                                            # dyn/cm**2=>atm) / Boltzmann's constant (R/avog)
        self._avog = 6.023e23 # Avogadro's number

        '''Setting simulation parameters'''
        self.iRep = 2   
        self.iStep = 0  # counts the number of vaporization steps
        self.iPrn = 0   # used to determine which steps are printed to the output file
        self.iIt = 0    # counts the number of iterations needed to solve for the activities
        self.iFirst = 0
        self.T = 2000   # temperature of magma (Kelvin)
        self.perVap = 0 # number of iterations of vaporization process

        '''Importing molecular oxide and metal weights'''
        mwOxides_fname = 'data/weights_oxides.csv'
        wMetals_fname = 'data/weights_metals.csv'

        self._mwOxides = dict(np.genfromtxt(mwOxides_fname,delimiter= ',',names=True,skip_header=1,\
                                 dtype=None, encoding = 'UTF-8'))
        self._wMetals = dict(np.genfromtxt(wMetals_fname,delimiter= ',',names=True,skip_header=1,\
                                 dtype=None, encoding = 'UTF-8',autostrip=True))

        # Evaluating the multiplications and dividing by _avog (mol wt. in g/mole)
        self._wMetals.update((x, eval(y)/self._avog) for x, y in self._wMetals.items()) 

    # end __init__()

    def start(self,model_comp):

        ''' Calculating number of moles of the oxides calculated from the input composition '''
        
        # Total weight
        self.totWt = sum(model_comp.wt_oxides.values())
        
        # Number of moles per oxide
        self.numMolOx = {}
        for ox in self._mwOxides:
            if model_comp.wt_oxides[ox] != 0: # Avoiding div by 0 error
                self.numMolOx[ox] = self._mwOxides[ox] / model_comp.wt_oxides[ox]
            else: 
                self.numMolOx[ox] = 0

        # Total amount of moles
        self.totMol = sum(self.numMolOx.values())

        # Mole percentages of the oxides
        self.mPerc = {}
        for ox in self.numMolOx:
            if model_comp.wt_oxides[ox] != 0: # Avoiding div by 0 error
                self.mPerc[ox] = self.numMolOx[ox] / self.totMol
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

    # end start()

    def set_temp(self,T):
        self.T = T
        print(f'Temperature set to {self.T} (K)')

    # end set_temp()

    def set_perVap(self,perVap):
        self.perVap = perVap
        print(f'Number of vaporization permutations set to {self.perVap}')

    # end set_perVap()

# end simulations()

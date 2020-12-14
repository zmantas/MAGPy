import numpy as np

class melt_activity():
    '''
    Calculates the activity of the oxides in the melt for a given composition
    and temperature using IMCC.
    '''
    def __init__(self,sim):
        
        '''
        The equilibrium constants for each of the relevant oxide pseudospecies
        is calculated below for the given temperature. The log10K values from
        which the A and B values are derived are given in the comments. See 
        Fegley et al. (1987) and Shaefer and Fegley (2004) for the source of 
        these values. 
        '''

        self.name_pseudo = {} # Names of oxide pseudospecies
        self.K_pseudo = {} # Equilibrium constant for oxide pseudospecies
        self.act_pseudo = {}
        '''
        MgO(liq) + SiO2(liq) = MgSiO3(liq)
        log10K(MgSiO3) = - 23.67 + 102856/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.name_pseudo['MG1'] = 'MgSiO3'
        self.K_pseudo['MG1'] = 10**( 0.42 + 2329/sim.T)

        '''
        2MgO(liq) + SiO2(liq) = Mg2SiO4(liq)
        log10K(Mg2SiO4) = - 34.08 + 141582/T
        -2log10K(MgO)   =  18.10 - 67242/T
         -log10K(SiO2)  =  15.04 - 66906/T
        '''
        self.name_pseudo['MG2'] = 'Mg2SiO4'
        self.K_pseudo['MG2'] = 10**(- 0.94 + 7434/sim.T)

        '''
        MgO(liq) + Al2O3(liq) = MgAl2O4(liq)
        log10K(MgAl2O4) = - 31.55 + 142219/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.name_pseudo['MG3'] = 'MgAl2O4'
        self.K_pseudo['MG3'] = 10**(1.18 + 464/sim.T)

        '''
        MgO(liq) + TiO2(liq) = MgTiO3(liq)
        log10K(MgTiO3) = - 22.54 + 103180/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.name_pseudo['MG4'] = 'MgTiO3'
        self.K_pseudo['MG4'] = 10**(- 0.13 + 3246/sim.T)

        '''
        MgO(liq) + 2TiO2(liq) = MgTi2O5(liq)
        log10K(MgTi2O5) = - 35.26 + 169092/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10(TiO2) = 26.72 - 132626/T
        '''
        self.name_pseudo['MG5'] = 'MgTi2O5'
        self.K_pseudo['MG5'] = 10**(0.51 + 2845/sim.T)

        '''
        2MgO(liq) + TiO2(liq) = Mg2TiO4(liq)
        log10K(Mg2TiO4) = - 30.79 + 137367/T
        -2log10K(MgO) = 18.10 - 67242/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.name_pseudo['MG6'] = 'Mg2TiO4'
        self.K_pseudo['MG6'] = 10**(0.67 + 3812/sim.T)

        '''
        2MgO(liq) + 2Al2O3(liq) + 5SiO2(liq) = Mg2Al4Si5O18(liq)
        log10K(Mg2Al4Si5O18) = - 132.38 + 618040/T
        -2log10K(MgO) = 18.10 - 67242/T
        -2log10K(Al2O3) = 47.36 - 216268/T
        -5log10K(SiO2) = 75.20 - 334530/T
        '''
        self.name_pseudo['MG7'] = 'Mg2Al4Si5O18'
        self.K_pseudo['MG7'] = 10**7.48

        '''
        3Al2O3(liq) + 2SiO2(liq) = Al6Si2O13(liq)
        log10K(Al6Si2O13) = - 104.06 + 467589/T
        -3log10K(Al2O3) = 71.04 - 324402/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.name_pseudo['AL1'] = 'Al6Si2O13'
        self.K_pseudo['AL1'] = 10**(- 2.94 + 9375/sim.T)

        '''
        CaO(liq) + Al2O3(liq) = CaAl2O4(liq)
        log10K(CaAl2O4) = - 33.93 + 154384/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.name_pseudo['CA1'] = 'CaAl2O4'
        self.K_pseudo['CA1'] = 10**(- 1.89 + 10060/sim.T)

        '''
        CaO(liq) + 2Al2O3(liq) = CaAl4O7(liq)
        log10K(CaAl4O7) = - 56.31 + 262171/T
        -log10K(CaO) = 8.36 - 36190/T
        -2log10K(Al2O3) = 47.36 - 216268/T
        '''
        self.name_pseudo['CA2'] = 'CaAl4O7'
        self.K_pseudo['CA2'] = 10**(- 0.59 + 9713/sim.T)

        '''
        12CaO(liq) + 7Al2O3(liq) = Ca12Al14O33(liq)
        log10K(Ca12Al14O33) = - 272.38 + 1263457/T
        -12log10K(CaO) = 100.32 - 434280/T
        -7log10K(Al2O3) = 165.76 - 756938/T
        '''
        self.name_pseudo['CA3'] = 'Ca12Al14O33'
        self.K_pseudo['CA3'] = 10**(-6.30 + 72239/sim.T)

        '''
        CaO(liq) + SiO2(liq) = CaSiO3(liq)
        log10K(CaSiO3) = - 22.86 + 108664/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.name_pseudo['CA4'] = 'CaSiO3'
        self.K_pseudo['CA4'] = 10**(0.54 + 5568/sim.T)

        '''
        CaO(liq) + Al2O3(liq) + 2SiO2(liq) = CaAl2Si2O8(liq)
        log10K(CaAl2Si2O8) = - 59.49 + 283462/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(Al2O3) = 23.68 - 108134/T
        -2log10(SiO2) = 30.08 - 133812/T
        '''
        self.name_pseudo['CA5'] = 'CaAl2Si2O8'
        self.K_pseudo['CA5'] = 10**(2.63 + 5326/sim.T)

        '''
        CaO(liq) + MgO(liq) + 2SiO2(liq) = CaMgSi2O6(liq)
        log10K(CaMgSi2O6) = -46.03 + 212108/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.name_pseudo['CA6'] = 'CaMgSi2O6'
        self.K_pseudo['CA6'] = 10**(1.46 + 8485/sim.T)

        '''
        2CaO(liq) + MgO(liq) + 2SiO2(liq) = Ca2MgSi2O7(liq)
        log10K(Ca2MgSi2O7) = - 55.22 + 255140/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.name_pseudo['CA7'] = 'Ca2MgSi2O7'
        self.K_pseudo['CA7'] = 10**(0.63 + 15327/sim.T)

        '''
        2CaO(liq) + Al2O3(liq) + SiO2(liq) = Ca2Al2SiO7(liq)
        log10K(Ca2Al2SiO7) = - 53.43 + 258130/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(Al2O3) = 23.68 - 108134/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.name_pseudo['CA8'] = 'Ca2Al2SiO7'
        self.K_pseudo['CA8'] = 10**(2.01 + 10710/sim.T)

        '''
        CaO(liq) + TiO2(liq) = CaTiO3(liq)
        log10K(CaTiO3) = - 21.80 + 109558/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.name_pseudo['CA9'] = 'CaTiO3'
        self.K_pseudo['CA9'] = 10**(- 0.08 + 7055/sim.T)

        '''
        2CaO(liq) + SiO2(liq) = Ca2SiO4(liq)
        log10K(Ca2SiO4) = - 31.13 + 147702/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.name_pseudo['CA10'] = 'Ca2SiO4'
        self.K_pseudo['CA10'] = 10**(0.63 + 8416/sim.T)

        '''
        CaO(liq) + TiO2(liq) + SiO2(liq) = CaTiSiO5(liq)
        log10K(CaTiSiO5) = - 36.94 + 179480/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(TiO2) = 13.36 - 66313/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.name_pseudo['CA11'] = 'CaTiSiO5'
        self.K_pseudo['CA11'] = 10**(- 0.18 + 10071/sim.T)

        '''
        CaO(liq) + 6Al2O3(liq) = CaAl12O19(liq)
        log10K(CaAl12O19) = -154.23 + 707606/T
        -log10K(CaO) = 8.36 - 36190/T
        -6log10K(Al2O3) = 142.08 - 648804/T
        '''
        self.name_pseudo['CA12'] = 'CaAl12O19'
        self.K_pseudo['CA12'] = 10**(- 3.79 + 22612/sim.T)

        '''
        FeO(liq) + TiO2(liq) = FeTiO3(liq)
        log10K(FeTiO3) = - 22.14 + 100392/T
        -log10K(FeO) = 8.27 - 30510/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.name_pseudo['FE1'] = 'FeTiO3'
        self.K_pseudo['FE1'] = 10**(- 0.51 + 3569/sim.T)
        # print('YOO',self.K_pseudo['FE1'])

        '''
        2FeO(liq) + SiO2(liq) = Fe2SiO4(liq)
        log10K(Fe2SiO4) = - 32.21 + 131029/T
        -2log10K(FeO) = 16.54 - 61020/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.name_pseudo['FE2'] = 'Fe2SiO4'
        self.K_pseudo['FE2'] = 10**(- 0.63 + 3103/sim.T)

        '''
        FeO(liq) + Al2O3(liq) = FeAl2O4(liq)
        log10K(FeAl2O4) = - 33.71 + 144336/T
        -log10K(FeO) = 8.27 - 30510/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.name_pseudo['FE3'] = 'FeAl2O4'
        self.K_pseudo['FE3'] = 10**(- 1.76 + 5692/sim.T)

        '''
        FeO (liq) + Fe2O3 (liq) = Fe3O4 (liq)
        Fe3O4 data from Barin 1995
        '''
        self.name_pseudo['FE4'] = 'Fe3O4'
        self.K_pseudo['FE4'] = 10**(-4.385894544e-1 + 4.3038155175436e3/sim.T\
                                   -3.1050205223386055e6/sim.T**2.0)

        '''
        Na2O(liq) + SiO2(liq) = Na2SiO3(liq)
        '''
        self.name_pseudo['NA1'] = 'Na2SiO3'
        self.K_pseudo['NA1'] = 10**(- 1.33 + 13870/sim.T)

        '''
        Na2O(liq) + 2SiO2(liq) = Na2Si2O5(liq)
        '''
        self.name_pseudo['NA2'] = 'Na2Si2O5'
        self.K_pseudo['NA2'] = 10**(- 1.39 + 15350/sim.T)

        ''' 
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + SiO2(liq) = NaAlSiO4(liq)
        '''
        self.name_pseudo['NA3'] = 'NaAlSiO4'
        self.K_pseudo['NA3'] = 10**(0.65 + 6997/sim.T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + 3SiO2(liq) = NaAlSi3O8(liq)
        '''
        self.name_pseudo['NA4'] = 'NaAlSi3O8'
        self.K_pseudo['NA4'] = 10**(1.29 + 8788/sim.T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) = NaAlO2(liq)
        '''
        self.name_pseudo['NA5'] = 'NaAlO2'
        self.K_pseudo['NA5'] = 10**(0.55 + 3058/sim.T)

        '''
        Na2O(liq) + TiO2(liq) = Na2TiO3(liq)
        '''
        self.name_pseudo['NA6'] = 'Na2TiO3'
        self.K_pseudo['NA6'] = 10**(- 1.38 + 15445/sim.T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + 2SiO2(liq) = NAAlSi2O6(liq)
        '''
        self.name_pseudo['NA7'] = 'NaAlSi2O6'
        self.K_pseudo['NA7'] = 10**(- 1.02 + 9607/sim.T)

        '''
        K2O(liq) + SiO2(liq) = K2SiO3(liq)
        '''
        self.name_pseudo['K1'] = 'K2SiO3'
        self.K_pseudo['K1'] = 10**(0.2692 + 12735/sim.T)

        '''
        K2O(liq) + 2SiO2(liq) = K2Si2O5(liq)
        '''
        self.name_pseudo['K2'] = 'K2Si2O5'
        self.K_pseudo['K2'] = 10**(0.3462 + 14685/sim.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + SiO2(liq) = KAlSiO4(liq)
        '''
        self.name_pseudo['K3'] = 'KAlSiO4'
        self.K_pseudo['K3'] = 10**(0.97 + 8675/sim.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + 3SiO2(liq) = KAlSi3O8(liq)
        '''
        self.name_pseudo['K4'] = 'KAlSi3O8'
        self.K_pseudo['K4'] = 10**(1.11 + 11229/sim.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) = KAlO2(liq)
        '''
        self.name_pseudo['K5'] = 'KAlO2'
        self.K_pseudo['K5'] = 10**(0.72 + 4679/sim.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + 2SiO2(liq) = KAlSi2O6(liq)
        '''
        self.name_pseudo['K6'] = 'KAlSi2O6'
        self.K_pseudo['K6'] = 10**(1.53 + 10125/sim.T)

        '''
        K2O(liq) + 4SiO2 (liq) = K2Si4O9 (liq)
        '''
        self.name_pseudo['K7'] = 'K2Si4O9'
        self.K_pseudo['K7'] = 10**(-0.9648 + 17572/sim.T)

        '''
        0.5K2O(liq) + CaO(liq) + 0.5Al2O3(liq) + 2SiO2(liq)=KCaAlSi2O7(liq)
        '''
        self.name_pseudo['K8'] = 'KCaAlSi2O7'
        self.K_pseudo['K8'] = 10**(4.2983 + 17037/sim.T)

    def activities_melt_pseudo(self,actOx):
        '''
        Calculates the activities for the complex species in the melt 
        (see activities_melt for relevant equations)
        '''
        self.act_pseudo['MG1'] = self.K_pseudo['MG1'] * actOx['MgO'] * actOx['SiO2']

        self.act_pseudo['MG2'] = self.K_pseudo['MG2'] * actOx['MgO']**2 * actOx['SiO2']
        
        self.act_pseudo['MG3'] = self.K_pseudo['MG3'] * actOx['MgO'] * actOx['Al2O3']
        
        self.act_pseudo['MG4'] = self.K_pseudo['MG4'] * actOx['MgO'] * actOx['TiO2']
        
        self.act_pseudo['MG5'] = self.K_pseudo['MG5'] * actOx['MgO'] * actOx['TiO2']**2
        
        self.act_pseudo['MG6'] = self.K_pseudo['MG6'] * actOx['MgO']**2 * actOx['TiO2']
        
        self.act_pseudo['MG7'] = self.K_pseudo['MG7'] * actOx['MgO']**2 * actOx['Al2O3']**2 * actOx['SiO2']**5
        
        self.act_pseudo['AL1'] = self.K_pseudo['AL1'] * actOx['Al2O3']**3 * actOx['SiO2']**2
        
        self.act_pseudo['CA1'] = self.K_pseudo['CA1'] * actOx['CaO'] * actOx['Al2O3']
        
        self.act_pseudo['CA2'] = self.K_pseudo['CA2'] * actOx['CaO'] * actOx['Al2O3']**2
        
        self.act_pseudo['CA3'] = self.K_pseudo['CA3'] * actOx['CaO']**12 * actOx['Al2O3']**7
        
        self.act_pseudo['CA4'] = self.K_pseudo['CA4'] * actOx['CaO'] * actOx['SiO2']
        
        self.act_pseudo['CA5'] = self.K_pseudo['CA5'] * actOx['CaO'] * actOx['Al2O3'] * actOx['SiO2']**2
        
        self.act_pseudo['CA6'] = self.K_pseudo['CA6'] * actOx['CaO'] * actOx['MgO'] * actOx['SiO2']**2
        
        self.act_pseudo['CA7'] = self.K_pseudo['CA7'] * actOx['CaO']**2 * actOx['MgO'] * actOx['SiO2']**2
        
        self.act_pseudo['CA8'] = self.K_pseudo['CA8'] * actOx['CaO']**2 * actOx['Al2O3'] * actOx['SiO2']
        
        self.act_pseudo['CA9'] = self.K_pseudo['CA9'] * actOx['CaO'] * actOx['TiO2']
        
        self.act_pseudo['CA10'] = self.K_pseudo['CA10'] * actOx['CaO']**2 * actOx['SiO2']
        
        self.act_pseudo['CA11'] = self.K_pseudo['CA11'] * actOx['CaO'] * actOx['TiO2'] * actOx['SiO2']
        
        self.act_pseudo['CA12'] = self.K_pseudo['CA12'] * actOx['CaO'] * actOx['Al2O3']**6

        self.act_pseudo['FE1'] = self.K_pseudo['FE1'] * actOx['FeO'] * actOx['TiO2']
          
        self.act_pseudo['FE2'] = self.K_pseudo['FE2'] * actOx['FeO']**2 * actOx['SiO2']
        
        self.act_pseudo['FE3'] = self.K_pseudo['FE3'] * actOx['FeO'] * actOx['Al2O3']
        
        self.act_pseudo['FE4'] = self.K_pseudo['FE4'] * actOx['FeO'] * actOx['Fe2O3']
                
        self.act_pseudo['NA1'] = self.K_pseudo['NA1'] * actOx['Na2O'] * actOx['SiO2']
        
        self.act_pseudo['NA2'] = self.K_pseudo['NA2'] * actOx['Na2O'] * actOx['SiO2']**2
        
        self.act_pseudo['NA3'] = self.K_pseudo['NA3'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']
        
        self.act_pseudo['NA4'] = self.K_pseudo['NA4'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**3
        
        self.act_pseudo['NA5'] = self.K_pseudo['NA5'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3'])
        
        self.act_pseudo['NA6'] = self.K_pseudo['NA6'] * actOx['Na2O'] * actOx['TiO2']
        
        self.act_pseudo['NA7'] = self.K_pseudo['NA7'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**2
        
        self.act_pseudo['K1'] = self.K_pseudo['K1'] * actOx['K2O'] * actOx['SiO2']
        
        self.act_pseudo['K2'] = self.K_pseudo['K2'] * actOx['K2O'] * actOx['SiO2']**2
        
        self.act_pseudo['K3'] = self.K_pseudo['K3'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']
        
        self.act_pseudo['K4'] = self.K_pseudo['K4'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**3
        
        self.act_pseudo['K5'] = self.K_pseudo['K5'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3'])
        
        self.act_pseudo['K6'] = self.K_pseudo['K6'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**2
        
        self.act_pseudo['K7'] = self.K_pseudo['K7'] * actOx['K2O'] * (actOx['SiO2']**4)
        
        self.act_pseudo['K8'] = self.K_pseudo['K8'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['CaO'] * actOx['SiO2']**2


    def recompute_gamma(self,act_ox,addF2O3):
        ''' 
        Recompute the activity coefficients of the oxides (gamma['Element']) from the acitivites
        computed above. 

        gamma['Element'] = 

        Activity(pure oxide) / SUM activities of all complex melt species containing the oxide
        
        Variables:
            - actOx: is required to calculate gamma for all elements
            - gamma: original value being passed on 
            - comp: necessary for calculating gamma['Fe2O3']
            - fAbMolecule: necessary for calculating gamma['Fe2O3']
        '''
        gamma_new = {} 

        #### gamma SiO2 ####
        if act_ox['SiO2'] != 0:
            gamma_new['SiO2'] = act_ox['SiO2'] / (act_ox['SiO2'] + \
                                self.act_pseudo['MG1'] + self.act_pseudo['MG2'] + \
                                self.act_pseudo['CA4'] + self.act_pseudo['CA8'] + \
                                self.act_pseudo['CA10']+ self.act_pseudo['CA11']+ \
                                self.act_pseudo['FE2'] + self.act_pseudo['NA1'] + \
                                self.act_pseudo['NA3'] + self.act_pseudo['K1']  + \
                                self.act_pseudo['K3']  + \
                                2 * (self.act_pseudo['AL1'] + self.act_pseudo['CA5'] + \
                                     self.act_pseudo['CA6'] + self.act_pseudo['CA7'] + \
                                     self.act_pseudo['NA2'] + self.act_pseudo['K2']  + \
                                     self.act_pseudo['K6']  + self.act_pseudo['K7']  + \
                                     self.act_pseudo['K8']) + \
                                3 * (self.act_pseudo['NA4'] + self.act_pseudo['K4']) + \
                                4 * (self.act_pseudo['K7']) + \
                                5 * (self.act_pseudo['MG7'])
                                )
        else:
            gamma_new['SiO2'] = 0

        #### gamma MgO ####
        if act_ox['MgO'] != 0:
            gamma_new['MgO'] = act_ox['MgO'] / (act_ox['MgO'] + \
                               self.act_pseudo['MG1'] + self.act_pseudo['MG3'] + \
                               self.act_pseudo['MG4'] + self.act_pseudo['MG5'] + \
                               self.act_pseudo['CA6'] + self.act_pseudo['CA7'] + \
                               2 * (self.act_pseudo['MG2'] + self.act_pseudo['MG6'] + \
                                    self.act_pseudo['MG7'])
                                 )
        else:
            gamma_new['MgO'] = 0
        
        #### gamma FeO #### 
        if act_ox['FeO'] != 0:
            gamma_new['FeO'] = act_ox['FeO'] / (act_ox['FeO'] + \
                                         self.act_pseudo['FE1'] + self.act_pseudo['FE3'] +\
                                         2 * (self.act_pseudo['FE2'] + act_ox['Fe2O3']) + \
                                         3 * (self.act_pseudo['FE4'])
                                         )
        else:
            gamma_new['FeO'] = 0

        #### gamma Fe3 ####
        # gamma['Fe3'] is an adjustment factor, not a true activity coefficient because
        # the mole fraction of Fe2O3 in the melt is not known.
        if not addF2O3:
            gamma_new['Fe2O3'] = 1
        
        elif addF2O3 and act_ox['Fe2O3'] != 0 and act_ox['FeO'] != 0:
            gamma_new['Fe2O3'] = act_ox['Fe2O3'] / (act_ox['Fe2O3'] + self.act_pseudo['FE4'])
        
        else:
            gamma_new['Fe2O3'] = 0
        
        #### gamma Ca ####
        if act_ox['CaO'] != 0:
            gamma_new['CaO'] = act_ox['CaO'] / (act_ox['CaO'] + \
                                         self.act_pseudo['CA1'] + self.act_pseudo['CA2'] + \
                                         self.act_pseudo['CA4'] + self.act_pseudo['CA5'] + \
                                         self.act_pseudo['CA6'] + self.act_pseudo['CA9'] + \
                                         self.act_pseudo['CA11']+ self.act_pseudo['CA12']+ \
                                         self.act_pseudo['K8']  + \
                                         2 * (self.act_pseudo['CA7'] + self.act_pseudo['CA8'] + \
                                              self.act_pseudo['CA10']) + \
                                         12 * (self.act_pseudo['CA3'])
                                         )
        else:
            gamma_new['CaO'] = 0
        
        #### gamma Al ####
        if act_ox['Al2O3'] != 0:
            gamma_new['Al2O3'] = act_ox['Al2O3'] / (act_ox['Al2O3'] + \
                                            self.act_pseudo['MG3'] + self.act_pseudo['CA1'] + \
                                            self.act_pseudo['CA5'] + self.act_pseudo['CA8'] + \
                                            self.act_pseudo['FE3'] + \
                                            2 * (self.act_pseudo['CA2'] + self.act_pseudo['MG7']) + \
                                            3 * (self.act_pseudo['AL1'])  + \
                                            6 * (self.act_pseudo['CA12']) + \
                                            7 * (self.act_pseudo['CA3'])  + \
                                            0.5 * (self.act_pseudo['NA3'] + self.act_pseudo['NA4'] + \
                                                   self.act_pseudo['NA5'] + self.act_pseudo['NA7'] + \
                                                   self.act_pseudo['K3']  + self.act_pseudo['K4']  + \
                                                   self.act_pseudo['K5']  + self.act_pseudo['K6']  + \
                                                   self.act_pseudo['K8'])
                                            )
        else: 
            gamma_new['Al2O3'] = 0
        
        #### gamma Ti ####
        if act_ox['TiO2'] != 0:
            gamma_new['TiO2'] = act_ox['TiO2'] / (act_ox['TiO2'] + \
                                           self.act_pseudo['MG4'] + self.act_pseudo['MG6'] + \
                                           self.act_pseudo['CA9'] + self.act_pseudo['CA11']+ \
                                           self.act_pseudo['FE1'] + self.act_pseudo['NA6'] + \
                                           2 * (self.act_pseudo['MG5'])
                                           )
        else:
            gamma_new['TiO2'] = 0
        
        #### gamma Na2 ####
        if act_ox['Na2O'] != 0:
            gamma_new['Na2O'] = act_ox['Na2O'] / (act_ox['Na2O'] + \
                                          self.act_pseudo['NA1'] + self.act_pseudo['NA2'] + \
                                          self.act_pseudo['NA6'] + \
                                          0.5 * (self.act_pseudo['NA3'] + self.act_pseudo['NA4'] + \
                                                 self.act_pseudo['NA5'] + self.act_pseudo['NA7'])
                                          )
        else:
            gamma_new['Na2O'] = 0
        
        #### gamma ###
        if act_ox['K2O'] != 0:
            gamma_new['K2O'] = act_ox['K2O'] / (act_ox['K2O'] + \
                                        self.act_pseudo['K1'] + self.act_pseudo['K2'] + \
                                        self.act_pseudo['K7'] + \
                                        0.5 * (self.act_pseudo['K8'] + self.act_pseudo['K3'] + \
                                               self.act_pseudo['K4'] + self.act_pseudo['K5'] + \
                                               self.act_pseudo['K6'])
                                        )
        else:
            gamma_new['K2O'] = 0

        return gamma_new

    def melt_activity_calculation(self,sim,addF2O3=False):

        iit = 0
        while iit < 1e8: 

            for metal in sim._metalNames:
                if metal != 'Fe3':
                    sim.act_ox[sim._metal2oxide[metal]] = sim.fAbOx[metal]\
                                        * sim.gamma[sim._metal2oxide[metal]]
                else:
                # Activity of Fe2O3 is estimated using gas chemistry, 
                # then all acitivities are recomputed
                    if addF2O3:
                        sim.act_ox['Fe2O3'] = sim.presLiq['Fe2O3']\
                                              * sim.gamma['Fe2O3']
                    else:
                        sim.act_ox['Fe2O3'] = 0
            
            # Compute activities of pseudo species 
            self.activities_melt_pseudo(sim.act_ox)

            # Recompute gamma with the updated activity values for the oxides
            gamma_new = self.recompute_gamma(sim.act_ox,addF2O3)

            # Compute ratio of newly computed activity and previous activity
            fGam = {}
            for oxide in sim.gamma: 
                if gamma_new[oxide] == 0 or \
                   (oxide == 'Fe2O3' and gamma_new['FeO'] == 0):
                    fGam[oxide] = 1
                else:
                    fGam[oxide] = gamma_new[oxide]/sim.gamma[oxide]

            # Could also be written in list comprenhension format, check numba:
            # self.fGam = {oxide : 1 if gamma_new[oxide] == 0 or\
            #            (oxide == 'Fe2O3' and gamma_new['FeO'] == 0) else\
            #            gamma_new[oxide]/sim.gamma[oxide]\
            #            for oxide in self.gama}
            
            ''' 
            If self.gamRat[elements] ~1, the code has arrived at a solution for all the 
            activities, and moves on to the gas chemistry. 
            If this is not the case, then the activity coefficients are adjusted and 
            the activities are recomputed until a solution is foudn. 
            ''' 

            if all(rat < 1e-5 for rat in np.abs(np.log10(list(fGam.values())))):
                break


            for element in gamma_new:
                if iit > 500:
                    sim.gamma[element] = (gamma_new[element] * sim.gamma[element]**4)**(1/5)
                elif iit > 30:
                    sim.gamma[element] = (gamma_new[element] * sim.gamma[element]**2)**(1/3)
                else:
                    sim.gamma[element] = (gamma_new[element] * sim.gamma[element])**(1/2)
            
            iit += 1 # updating counter
            
            if iit >= 1e8: 
                raise RuntimeError('Max recursion limit reached while calculating activities.')


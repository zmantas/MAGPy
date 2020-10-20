'''
Gas chemistry thermodynamic data
'''

import numpy as np

class thermodynamic_data():
   
    def __init__(self,T):
       
        self.T = T   # temperature
        self.actMelt = {} # melt activities
        self.actMeltComp = {}
        self.nameMeltComp = {}

    # end __init__()

    def activities_melt(self):
        '''
        Thermodynamic data for activities in the melt
            returns melt activities for given temperature
        '''
        '''
        MgO(liq) + SiO2(liq) = MgSiO3(liq)
        log10K(MgSiO3) = - 23.67 + 102856/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.nameMeltComp['MG1'] = 'MgSiO3'
        self.actMelt['MG1'] = 10**( 0.42 + 2329/self.T)

        '''
        2MgO(liq) + SiO2(liq) = Mg2SiO4(liq)
        log10K(Mg2SiO4) = - 34.08 + 141582/T
        -2log10K(MgO)   =  18.10 - 67242/T
         -log10K(SiO2)  =  15.04 - 66906/T
        '''
        self.nameMeltComp['MG2'] = 'Mg2SiO4'
        self.actMelt['MG2'] = 10**(- 0.94 + 7434/self.T)

        '''
        MgO(liq) + Al2O3(liq) = MgAl2O4(liq)
        log10K(MgAl2O4) = - 31.55 + 142219/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.nameMeltComp['MG3'] = 'MgAl2O4'
        self.actMelt['MG3'] = 10**(1.18 + 464/self.T)

        '''
        MgO(liq) + TiO2(liq) = MgTiO3(liq)
        log10K(MgTiO3) = - 22.54 + 103180/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.nameMeltComp['MG4'] = 'MgTiO3'
        self.actMelt['MG4'] = 10**(- 0.13 + 3246/self.T)

        '''
        MgO(liq) + 2TiO2(liq) = MgTi2O5(liq)
        log10K(MgTi2O5) = - 35.26 + 169092/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10(TiO2) = 26.72 - 132626/T
        '''
        self.nameMeltComp['MG5'] = 'MgTi2O5'
        self.actMelt['MG5'] = 10**(0.51 + 2845/self.T)

        '''
        2MgO(liq) + TiO2(liq) = Mg2TiO4(liq)
        log10K(Mg2TiO4) = - 30.79 + 137367/T
        -2log10K(MgO) = 18.10 - 67242/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.nameMeltComp['MG6'] = 'Mg2TiO4'
        self.actMelt['MG6'] = 10**(0.67 + 3812/self.T)

        '''
        2MgO(liq) + 2Al2O3(liq) + 5SiO2(liq) = Mg2Al4Si5O18(liq)
        log10K(Mg2Al4Si5O18) = - 132.38 + 618040/T
        -2log10K(MgO) = 18.10 - 67242/T
        -2log10K(Al2O3) = 47.36 - 216268/T
        -5log10K(SiO2) = 75.20 - 334530/T
        '''
        self.nameMeltComp['MG7'] = 'Mg2Al4Si5O18'
        self.actMelt['MG7'] = 10**7.48

        '''
        3Al2O3(liq) + 2SiO2(liq) = Al6Si2O13(liq)
        log10K(Al6Si2O13) = - 104.06 + 467589/T
        -3log10K(Al2O3) = 71.04 - 324402/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.nameMeltComp['AL1'] = 'Al6Si2O13'
        self.actMelt['AL1'] = 10**(- 2.94 + 9375/self.T)

        '''
        CaO(liq) + Al2O3(liq) = CaAl2O4(liq)
        log10K(CaAl2O4) = - 33.93 + 154384/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.nameMeltComp['CA1'] = 'CaAl2O4'
        self.actMelt['CA1'] = 10**(- 1.89 + 10060/self.T)

        '''
        CaO(liq) + 2Al2O3(liq) = CaAl4O7(liq)
        log10K(CaAl4O7) = - 56.31 + 262171/T
        -log10K(CaO) = 8.36 - 36190/T
        -2log10K(Al2O3) = 47.36 - 216268/T
        '''
        self.nameMeltComp['CA2'] = 'CaAl4O7'
        self.actMelt['CA2'] = 10**(- 0.59 + 9713/self.T)

        '''
        12CaO(liq) + 7Al2O3(liq) = Ca12Al14O33(liq)
        log10K(Ca12Al14O33) = - 272.38 + 1263457/T
        -12log10K(CaO) = 100.32 - 434280/T
        -7log10K(Al2O3) = 165.76 - 756938/T
        '''
        self.nameMeltComp['CA3'] = 'Ca12Al14O33'
        self.actMelt['CA3'] = 10**(-6.30 + 72239/self.T)

        '''
        CaO(liq) + SiO2(liq) = CaSiO3(liq)
        log10K(CaSiO3) = - 22.86 + 108664/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.nameMeltComp['CA4'] = 'CaSiO3'
        self.actMelt['CA4'] = 10**(0.54 + 5568/self.T)

        '''
        CaO(liq) + Al2O3(liq) + 2SiO2(liq) = CaAl2Si2O8(liq)
        log10K(CaAl2Si2O8) = - 59.49 + 283462/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(Al2O3) = 23.68 - 108134/T
        -2log10(SiO2) = 30.08 - 133812/T
        '''
        self.nameMeltComp['CA5'] = 'CaAl2Si2O8'
        self.actMelt['CA5'] = 10**(+ 2.63 + 5326/self.T)

        '''
        CaO(liq) + MgO(liq) + 2SiO2(liq) = CaMgSi2O6(liq)
        log10K(CaMgSi2O6) = -46.03 + 212108/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.nameMeltComp['CA6'] = 'CaMgSi2O6'
        self.actMelt['CA6'] = 10**(1.46 + 8485/self.T)

        '''
        2CaO(liq) + MgO(liq) + 2SiO2(liq) = Ca2MgSi2O7(liq)
        log10K(Ca2MgSi2O7) = - 55.22 + 255140/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.nameMeltComp['CA7'] = 'Ca2MgSi2O7'
        self.actMelt['CA7'] = 10**(0.63 + 15327/self.T)

        '''
        2CaO(liq) + Al2O3(liq) + SiO2(liq) = Ca2Al2SiO7(liq)
        log10K(Ca2Al2SiO7) = - 53.43 + 258130/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(Al2O3) = 23.68 - 108134/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.nameMeltComp['CA8'] = 'Ca2Al2SiO7'
        self.actMelt['CA8'] = 10**(2.01 + 10710/self.T)

        '''
        CaO(liq) + TiO2(liq) = CaTiO3(liq)
        log10K(CaTiO3) = - 21.80 + 109558/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.nameMeltComp['CA9'] = 'CaTiO3'
        self.actMelt['CA9'] = 10**(- 0.08 + 7055/self.T)

        '''
        2CaO(liq) + SiO2(liq) = Ca2SiO4(liq)
        log10K(Ca2SiO4) = - 31.13 + 147702/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.nameMeltComp['CA10'] = 'Ca2SiO4'
        self.actMelt['CA10'] = 10**(0.63 + 8416/self.T)

        '''
        CaO(liq) + TiO2(liq) + SiO2(liq) = CaTiSiO5(liq)
        log10K(CaTiSiO5) = - 36.94 + 179480/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(TiO2) = 13.36 - 66313/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.nameMeltComp['CA11'] = 'CaTiSiO5'
        self.actMelt['CA11'] = 10**(- 0.18 + 10071/self.T)

        '''
        CaO(liq) + 6Al2O3(liq) = CaAl12O19(liq)
        log10K(CaAl12O19) = -154.23 + 707606/T
        -log10K(CaO) = 8.36 - 36190/T
        -6log10K(Al2O3) = 142.08 - 648804/T
        '''
        self.nameMeltComp['CA12'] = 'CaAl12O19'
        self.actMelt['CA12'] = 10**(- 3.79 + 22612/self.T)

        '''
        FeO(liq) + TiO2(liq) = FeTiO3(liq)
        log10K(FeTiO3) = - 22.14 + 100392/T
        -log10K(FeO) = 8.27 - 30510/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.nameMeltComp['FE1'] = 'FeTiO3'
        self.actMelt['FE1'] = 10**(- 0.51 + 3569/self.T)

        '''
        2FeO(liq) + SiO2(liq) = Fe2SiO4(liq)
        log10K(Fe2SiO4) = - 32.21 + 131029/T
        -2log10K(FeO) = 16.54 - 61020/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.nameMeltComp['FE2'] = 'Fe2SiO4'
        self.actMelt['FE2'] = 10**(- 0.63 + 3103/self.T)

        '''
        FeO(liq) + Al2O3(liq) = FeAl2O4(liq)
        log10K(FeAl2O4) = - 33.71 + 144336/T
        -log10K(FeO) = 8.27 - 30510/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.nameMeltComp['FE3'] = 'FeAl2O4'
        self.actMelt['FE3'] = 10**(- 1.76 + 5692/self.T)

        '''
        FeO (liq) + Fe2O3 (liq) = Fe3O4 (liq)
        Fe3O4 data from Barin 1995
        '''
        self.nameMeltComp['FE4'] = 'Fe3O4'
        self.actMelt['FE4'] = 10**(-4.385894544 + 4.3038155175436/self.T - 3.1050205223386055/self.T**2.0)

        '''
        Na2O(liq) + SiO2(liq) = Na2SiO3(liq)
        '''
        self.nameMeltComp['NA1'] = 'Na2SiO3'
        self.actMelt['NA1'] = 10**(- 1.33 + 13870/self.T)

        '''
        Na2O(liq) + 2SiO2(liq) = Na2Si2O5(liq)
        '''
        self.nameMeltComp['NA2'] = 'Na2Si2O5'
        self.actMelt['NA2'] = 10**(- 1.39 + 15350/self.T)

        ''' 
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + SiO2(liq) = NaAlSiO4(liq)
        '''
        self.nameMeltComp['NA3'] = 'NaAlSiO4'
        self.actMelt['NA3'] = 10**(0.65 + 6997/self.T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + 3SiO2(liq) = NaAlSi3O8(liq)
        '''
        self.nameMeltComp['NA4'] = 'NaAlSi3O8'
        self.actMelt['NA4'] = 10**(1.29 + 8788/self.T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) = NaAlO2(liq)
        '''
        self.nameMeltComp['NA5'] = 'NaAlO2'
        self.actMelt['NA5'] = 10**(0.55 + 3058/self.T)

        '''
        Na2O(liq) + TiO2(liq) = Na2TiO3(liq)
        '''
        self.nameMeltComp['NA6'] = 'Na2TiO3'
        self.actMelt['NA6'] = 10**(- 1.38 + 15445/self.T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + 2SiO2(liq) = NAAlSi2O6(liq)
        '''
        self.nameMeltComp['NA7'] = 'NaAlSi2O6'
        self.actMelt['NA7'] = 10**(- 1.02 + 9607/self.T)

        '''
        K2O(liq) + SiO2(liq) = K2SiO3(liq)
        '''
        self.nameMeltComp['K1'] = 'K2SiO3'
        self.actMelt['K1'] = 10**(0.2692 + 12735/self.T)

        '''
        K2O(liq) + 2SiO2(liq) = K2Si2O5(liq)
        '''
        self.nameMeltComp['K2'] = 'K2Si2O5'
        self.actMelt['K2'] = 10**(0.3462 + 14685/self.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + SiO2(liq) = KAlSiO4(liq)
        '''
        self.nameMeltComp['K3'] = 'KAlSiO4'
        self.actMelt['K3'] = 10**(0.97 + 8675/self.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + 3SiO2(liq) = KAlSi3O8(liq)
        '''
        self.nameMeltComp['K4'] = 'KAlSi3O8'
        self.actMelt['K4'] = 10**(1.11 + 11229/self.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) = KAlO2(liq)
        '''
        self.nameMeltComp['K5'] = 'KAlO2'
        self.actMelt['K5'] = 10**(0.72 + 4679/self.T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + 2SiO2(liq) = KAlSi2O6(liq)
        '''
        self.nameMeltComp['K6'] = 'KAlSi2O6'
        self.actMelt['K6'] = 10**(1.53 + 10125/self.T)

        '''
        K2O(liq) + 4SiO2 (liq) = K2Si4O9 (liq)
        '''
        self.nameMeltComp['K7'] = 'K2Si4O9'
        self.actMelt['K7'] = 10**(-0.9648 + 17572/self.T)

        '''
        0.5K2O(liq) + CaO(liq) + 0.5Al2O3(liq) + 2SiO2(liq)=KCaAlSi2O7(liq)
        '''
        self.nameMeltComp['K8'] = 'KCaAlSi2O7'
        self.actMelt['K8'] = 10**(4.2983 + 17037/self.T)

        return self.actMelt

    # end activities_melt()

    def activities_meltComplex(self,actOx):
        '''
        Calculates the activities for the complex species in the melt 
        (see activities_melt for relevant equations)
        '''
        self.actMeltComp['MG1'] = self.actMelt['MG1'] * actOx['MgO'] * actOx['SiO2']
        
        self.actMeltComp['MG2'] = self.actMelt['MG2'] * actOx['MgO']**2 * actOx['SiO2']
        
        self.actMeltComp['MG3'] = self.actMelt['MG3'] * actOx['MgO'] * actOx['Al2O3']
        
        self.actMeltComp['MG4'] = self.actMelt['MG4'] * actOx['MgO'] * actOx['TiO2']
        
        self.actMeltComp['MG5'] = self.actMelt['MG5'] * actOx['MgO'] * actOx['TiO2']**2
        
        self.actMeltComp['MG6'] = self.actMelt['MG6'] * actOx['MgO']**2 * actOx['TiO2']
        
        self.actMeltComp['MG7'] = self.actMelt['MG7'] * actOx['MgO']**2 * actOx['Al2O3']**2 * actOx['SiO2']**5
        
        self.actMeltComp['AL1'] = self.actMelt['AL1'] * actOx['Al2O3']**3 * actOx['SiO2']**2
        
        self.actMeltComp['CA1'] = self.actMelt['CA1'] * actOx['CaO'] * actOx['Al2O3']
        
        self.actMeltComp['CA2'] = self.actMelt['CA2'] * actOx['CaO'] * actOx['Al2O3']**2
        
        self.actMeltComp['CA3'] = self.actMelt['CA3'] * actOx['CaO']**12 * actOx['Al2O3']**7
        
        self.actMeltComp['CA4'] = self.actMelt['CA4'] * actOx['CaO'] * actOx['SiO2']
        
        self.actMeltComp['CA5'] = self.actMelt['CA5'] * actOx['CaO'] * actOx['Al2O3'] * actOx['SiO2']**2
        
        self.actMeltComp['CA6'] = self.actMelt['CA6'] * actOx['CaO'] * actOx['MgO'] * actOx['SiO2']**2
        
        self.actMeltComp['CA7'] = self.actMelt['CA7'] * actOx['CaO']**2 * actOx['MgO'] * actOx['SiO2']**2
        
        self.actMeltComp['CA8'] = self.actMelt['CA8'] * actOx['CaO']**2 * actOx['Al2O3'] * actOx['SiO2']
        
        self.actMeltComp['CA9'] = self.actMelt['CA9'] * actOx['CaO'] * actOx['TiO2']
        
        self.actMeltComp['CA10'] = self.actMelt['CA10'] * actOx['CaO']**2 * actOx['SiO2']
        
        self.actMeltComp['CA11'] = self.actMelt['CA11'] * actOx['CaO'] * actOx['TiO2'] * actOx['SiO2']
        
        self.actMeltComp['CA12'] = self.actMelt['CA12'] * actOx['CaO'] * actOx['Al2O3']**6

        self.actMeltComp['FE1'] = self.actMelt['FE1'] * actOx['FeO'] * actOx['TiO2']
        
        self.actMeltComp['FE2'] = self.actMelt['FE2'] * actOx['FeO']**2 * actOx['SiO2']
        
        self.actMeltComp['FE3'] = self.actMelt['FE3'] * actOx['FeO'] * actOx['Al2O3']
        
        self.actMeltComp['FE4'] = self.actMelt['FE4'] * actOx['FeO'] * actOx['Fe2O3']
                
        self.actMeltComp['NA1'] = self.actMelt['NA1'] * actOx['Na2O'] * actOx['SiO2']
        
        self.actMeltComp['NA2'] = self.actMelt['NA2'] * actOx['Na2O'] * actOx['SiO2']**2
        
        self.actMeltComp['NA3'] = self.actMelt['NA3'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']
        
        self.actMeltComp['NA4'] = self.actMelt['NA4'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**3
        
        self.actMeltComp['NA5'] = self.actMelt['NA5'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3'])
        
        self.actMeltComp['NA6'] = self.actMelt['NA6'] * actOx['Na2O'] * actOx['TiO2']
        
        self.actMeltComp['NA7'] = self.actMelt['NA7'] * np.sqrt(actOx['Na2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**2
        
        self.actMeltComp['K1'] = self.actMelt['K1'] * actOx['K2O'] * actOx['SiO2']
        
        self.actMeltComp['K2'] = self.actMelt['K2'] * actOx['K2O'] * actOx['SiO2']**2
        
        self.actMeltComp['K3'] = self.actMelt['K3'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']
        
        self.actMeltComp['K4'] = self.actMelt['K4'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**3
        
        self.actMeltComp['K5'] = self.actMelt['K5'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3'])
        
        self.actMeltComp['K6'] = self.actMelt['K6'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['SiO2']**2
        
        self.actMeltComp['K7'] = self.actMelt['K7'] * actOx['K2O'] * (actOx['SiO2']**4)
        
        self.actMeltComp['K8'] = self.actMelt['K8'] * np.sqrt(actOx['K2O']) * np.sqrt(actOx['Al2O3']) * actOx['CaO'] * actOx['SiO2']**2

        return self.actMeltComp

    # end activities_meltComplex()

    def recompute_gamma(self, actOx, gamma, comp, fAbMolecule):

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

        TODO: This part of the code gives slightly different values than the F77 code.
              This is most likely due to a difference in precesion. Should be checked out. 

        '''
        gamma = gamma.copy() # Necessary in order to avoid working in original gamma value dict

        #### gamma SiO2 ####
        if actOx['SiO2'] != 0:
            gamma['Si'] = actOx['SiO2'] / (actOx['SiO2'] + \
                                           self.actMeltComp['MG1'] + self.actMeltComp['MG2'] + \
                                           self.actMeltComp['CA4'] + self.actMeltComp['CA8'] + \
                                           self.actMeltComp['CA10']+ self.actMeltComp['CA11']+ \
                                           self.actMeltComp['FE2'] + self.actMeltComp['NA1'] + \
                                           self.actMeltComp['NA3'] + self.actMeltComp['K1']  + \
                                           self.actMeltComp['K3']  + \
                                           2 * (self.actMeltComp['AL1'] + self.actMeltComp['CA5'] + \
                                                self.actMeltComp['CA6'] + self.actMeltComp['CA7'] + \
                                                self.actMeltComp['NA2'] + self.actMeltComp['K2']  + \
                                                self.actMeltComp['K6']  + self.actMeltComp['K7']  + \
                                                self.actMeltComp['K8']) + \
                                           3 * (self.actMeltComp['NA4'] + self.actMeltComp['K4']) + \
                                           4 * (self.actMeltComp['K7']) + \
                                           5 * (self.actMeltComp['MG7'])
                                           )
        else:
            gamma['Si'] = 0

        #### gamma Mg ####
        if actOx['MgO'] != 0:
            gamma['Mg'] = actOx['MgO'] / (actOx['MgO'] + \
                                         self.actMeltComp['MG1'] + self.actMeltComp['MG3'] + \
                                         self.actMeltComp['MG4'] + self.actMeltComp['MG5'] + \
                                         self.actMeltComp['CA6'] + self.actMeltComp['CA7'] + \
                                         2 * (self.actMeltComp['MG2'] + self.actMeltComp['MG6'] + \
                                              self.actMeltComp['MG7'])
                                         )
        else:
            gamma['Mg'] = 0
        
        #### gamma Fe #### 
        if actOx['FeO'] != 0:
            gamma['Fe'] = actOx['FeO'] / (actOx['FeO'] + \
                                         self.actMeltComp['FE1'] + self.actMeltComp['FE3'] +\
                                         2 * (self.actMeltComp['FE2'] + actOx['Fe2O3']) + \
                                         3 * (self.actMeltComp['FE4'])
                                         )
        else:
            gamma['Fe'] = 0

        #### gamma Fe3 ####
        # gamma['Fe3'] is an adjustment factor, not a true activity coefficient because
        # the mole fraction of Fe2O3 in the melt is not known.
        if comp == 1:
            gamma['Fe3'] = 1
        
        elif comp != 1 and actOx['Fe2O3'] != 0 and fAbMolecule['Fe'] != 0:
            gamma['Fe3'] = actOx['Fe2O3'] / (actOx['Fe2O3'] + self.actMeltComp['FE4'])
        
        else:
            gamma['Fe3'] = 0
        
        #### gamma Ca ####
        if actOx['CaO'] != 0:
            gamma['Ca'] = actOx['CaO'] / (actOx['CaO'] + \
                                         self.actMeltComp['CA1'] + self.actMeltComp['CA2'] + \
                                         self.actMeltComp['CA4'] + self.actMeltComp['CA5'] + \
                                         self.actMeltComp['CA6'] + self.actMeltComp['CA9'] + \
                                         self.actMeltComp['CA11']+ self.actMeltComp['CA12']+ \
                                         self.actMeltComp['K8']  + \
                                         2 * (self.actMeltComp['CA7'] + self.actMeltComp['CA8'] + \
                                              self.actMeltComp['CA10']) + \
                                         12 * (self.actMeltComp['CA3'])
                                         )
        else:
            gamma['Ca'] = 0
        
        #### gamma Al ####
        if actOx['Al2O3'] != 0:
            gamma['Al'] = actOx['Al2O3'] / (actOx['Al2O3'] + \
                                            self.actMeltComp['MG3'] + self.actMeltComp['CA1'] + \
                                            self.actMeltComp['CA5'] + self.actMeltComp['CA8'] + \
                                            self.actMeltComp['FE3'] + \
                                            2 * (self.actMeltComp['CA2'] + self.actMeltComp['MG7']) + \
                                            3 * (self.actMeltComp['AL1'])  + \
                                            6 * (self.actMeltComp['CA12']) + \
                                            7 * (self.actMeltComp['CA3'])  + \
                                            0.5 * (self.actMeltComp['NA3'] + self.actMeltComp['NA4'] + \
                                                   self.actMeltComp['NA5'] + self.actMeltComp['NA7'] + \
                                                   self.actMeltComp['K3']  + self.actMeltComp['K4']  + \
                                                   self.actMeltComp['K5']  + self.actMeltComp['K6']  + \
                                                   self.actMeltComp['K8'])
                                            )
        else: 
            gamma['Al'] = 0
        
        #### gamma Ti ####
        if actOx['TiO2'] != 0:
            gamma['Ti'] = actOx['TiO2'] / (actOx['TiO2'] + \
                                           self.actMeltComp['MG4'] + self.actMeltComp['MG6'] + \
                                           self.actMeltComp['CA9'] + self.actMeltComp['CA11']+ \
                                           self.actMeltComp['FE1'] + self.actMeltComp['NA6'] + \
                                           2 * (self.actMeltComp['MG5'])
                                           )
        else:
            gamma['Ti'] = 0
        
        #### gamma Na2 ####
        if actOx['Na2O'] != 0:
            gamma['Na'] = actOx['Na2O'] / (actOx['Na2O'] + \
                                          self.actMeltComp['NA1'] + self.actMeltComp['NA2'] + \
                                          self.actMeltComp['NA6'] + \
                                          0.5 * (self.actMeltComp['NA3'] + self.actMeltComp['NA4'] + \
                                                 self.actMeltComp['NA5'] + self.actMeltComp['NA7'])
                                          )
        else:
            gamma['Na'] = 0
        
        #### gamma ###
        if actOx['K2O'] != 0:
            gamma['K'] = actOx['K2O'] / (actOx['K2O'] + \
                                        self.actMeltComp['K1'] + self.actMeltComp['K2'] + \
                                        self.actMeltComp['K7'] + \
                                        0.5 * (self.actMeltComp['K8'] + self.actMeltComp['K3'] + \
                                               self.actMeltComp['K4'] + self.actMeltComp['K5'] + \
                                               self.actMeltComp['K6'])
                                        )
        else:
            gamma['K'] = 0

        return gamma
        
    # end recompute_gamma()
        
# end melt_activities():
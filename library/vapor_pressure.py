import numpy as np

class vapor_pressure():

    def __init__(self,sim):

        ''' Setting initial values of key pressures and adjustment factors for gases '''
        # self.presGas = {gas : 1 for gas in sim._gasNames}  # gas pressures
        self.adjFact = {gas : 1 for gas in sim._gasNames}  # adjustment factors
        self._pConv = 1.01325e6/1.38046e-16 # _pConv converts the pressures into number         densities
                                            # dyn/cm**2=>atm) / Boltzmann's constant (R/avog)
        '''
        Gas chemistry thermodynamic data
        - This has been taken almost directly from the Fortran code due to lack of insight
          into the variable naming. Maybe derserves another look in the future.
        '''

        ''' 
        ###### Silicon and oxygen ######
        SiO2(liq) = Si(g) + 2O(g)
        JANAF 2nd ed. & supplements 2000-4500 K every 500 degrees 
        correlation coefficient for linear fit = -0.99989
        AK1 = 10.0**A = P(SIG)*P(OG)**2/P(SIO2L)
        '''    
        A = 22.13 - 94311.0 /sim.T
         
        '''
        0.5O2 (g) = O(g)
        JANAF 2nd ed. 
        correlation coefficient for linear fit = 
        AK2 = 10.0**E = P(OG)/P(02G)**0.5

        HENCE 10.0**A = P(SIG)*(P(O2G)*10.0**2E)/P(SIO2L)
        HENCE P(SIO2L) = 10.0**(2.0*E-A)*P(SIG)*P(O2G)
        HENCE P(O2G) = 10.0**(-2.0*E) * P(OG)**2
        '''
        E = 3.47 - 13282 /sim.T
        self.EOG = 10.0**E
        # print(self.EOG)
        
        '''    
        Si(liq) = Si(g)
        AK3 = 10.0**B = P(SIG)/P(SIL)
        HENCE P(SIG) = 10.0**B*P(SIL)
        '''
        B = 6.00 - 20919 /sim.T
        self.ESIG = 10.0**B

        '''
        Si(liq) + 0.5 O2(g) = SiO(g)
        FOR P(SIG) INSTEAD OF P(SIL)    C2 = - 3.67 + 29760 /sim.T
        AK4 = 10.0**C = P(SIOG)/(P(SIL)*P(O2G)**0.5)
        HENCE P(SIL) = 10.0**-C*P(SIOG)/P(O2G)**0.5
        HENCE P(SIO2L) = 10.0**(2.0*E+B-A)*10.0**(-C)*P(SIOG)*P(O2G)**0.5
        HENCE P(SIO2L) = 10.0**(2.0*E+B-A-C)*P(SIOG)*P(O2G)**0.5
        '''
        C = 2.51 + 8207 /sim.T
        self.ESIL = 10.0**(-C)
        self.ESIO2L = 10.0**(2.0*E + B - A - C)

        '''
        Si(liq) + O2(g) = SiO2(g)
        FOR P(SIG) INSTEAD OF P(SIL)    D2 = - 7.57 + 39595 /sim.T
        AK5 = 10.0**D = P(SIO2G)/(P(SIL)*P(O2G))
        HENCE P(SIO2G) = 10.0**D*P(SIL)*P(O2G)
        HENCE P(SIO2G) = 10.0**D*10.0**-B*P(SIG)*P(O2G)
        HENCE P(SIO2G) = 10.0**(D-B) * P(SIG) * P(O2G)
        '''
        D = -1.44 + 18326 /sim.T
        self.ESIO2G = 10.0**D

        '''
        ###### Magnesium ######
        MgO(liq) = Mg(g) + O(g)
        JANAF 2nd ed. & supplements 1500-5000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99994
        AK6 = 10.0**F = P(MGG)*P(OG)/P(MGOL)
        HENCE P(MGOL) = 10.0**(-F) * P(MGG) * P(OG)
        '''
        F = 12.56 - 46992 /sim.T
        self.EMGOL = 10.0**(-F)
    
        '''
        Mg(g) + 0.5 O2(g) = MgO(g)
        AK8 = 10.0**G = P(MGOG)/(P(MGG)*P(O2G)**0.5)
        HENCE P(MGOG) = 10.0**G * P(MGG) * P(O2G)**0.5
        '''
        G = -1.19 + 3794 /sim.T
        self.EMGG = 10.0**(-G)

        '''
        ###### Iron ######
        FeO(liq) = Fe(g) + O(g)
        JANAF 2nd ed. & supplements 2000-5000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99995
        AK9 = 10.0**AA = P(FEG)*P(OG)/P(FEOL)
        HENCE P(FEOL) = 10.0**(-AA) * P(FEG) * P(OG)
        '''
        AA = 12.06 - 44992 /sim.T
        self.EFEOL = 10.0**(-AA)
        
        '''
        Fe(liq) = Fe(g)
        AK10 = 10.0**AB = P(FEG)/P(FEL)
        Fe(liq) + 0.5 O2(g) = FeO(g)
        FOR P(FEG) INSTEAD OF P(FEL)    AC2 = - 2.93 + 9945 /sim.T
        AK11 = 10.0**AC = P(FEOG)/(P(FEL)*P(O2G)**0.5)
        HENCE P(FEOG) = 10.0**(AC-AB) * P(FEG) * P(O2G)**0.5
        '''
        AB = 6.35 -19704 /sim.T
        self.EFEL = 10.0**(-AB)
        AC = 3.39 - 9951 /sim.T
        self.EFEOG = 10.0**(AC-AB)

        '''
        2Fe (g) + 1.5 O2 (g) = Fe2O3 (liq)
        Fe2O3 (liq) cp data from IVTANTHERMO database (estimated)
        hematite enthalpy of fusion calculated from Sugawara & Akaogi 2004
        K = PFE2O3L / (PFEG**2 * PO2G**0.5)
        HENCE P(FE2O3L) = 10**AD * PFEG**2 * PO2G**0.5
        '''
        AD = -2.26722053113e1 + 7.56430936141329e4 /sim.T
        self.EFE2O3L = 10**AD
        
        '''
        3Fe (g) + 2 O2 (g) = Fe3O4 (liq)
        Fe3O4 (liq) cp data from Barin 95
        magnetite enthalpy of fusion from JANAF 4th ed.
        K = PFE3O4L / (PFEG**3 * PO2G**0.5)
        HENCE P(FE3O4L) = 10**AE * PFEG**3 * PO2G**0.5
        '''
        AE = -3.19907301154e1 + 1.110526206139634e5 /sim.T
        self.EFE3O4L = 10**AE

        '''
        ###### Calcium ######
        CaO(liq) = Ca(g) + O(g)
        JANAF 2nd ed. & supplements 2000-4500 K every 500 degrees 
        correlation coefficient for linear fit = -0.99998
        AK12 = 10.0**BA = P(CAG)*P(OG)/P(CAOL)
        HENCE P(CAOL) = 10.0**(-BA) * P(CAG) * P(OG)
        '''
        BA = 11.88 - 49586 /sim.T
        self.ECAOL = 10.0**(-BA)
        
        '''
        Ca(g) + 0.5 O2(g) = CaO(g)
        AK14 = 10.0**BC = P(CAOG)/(P(CAG)*P(O2G)**0.5)
        HENCE P(CAOG) = 10.0**BC * P(CAG) * P(O2G)**0.5
        '''
        BC = -1.61 + 6128 /sim.T
        self.ECAOG = 10.0**BC

        '''
        ###### Aluminum ###### 

        Al2O3(liq) = 2Al(g) + 3O(g)
        JANAF supplements 1500-4000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99995
        AK15 = 10.0**CA = P(ALG)**2*P(OG)**3/P(AL2O3L)
        HENCE P(AL2O3L) = 10.0**(-CA) * P(ALG)**2 * P(OG)**3
        '''
        CA = 35.83 - 153255 /sim.T
        self.EAL2O3L = 10.0**(-CA)
    
        '''
        Al(liq) = Al(g)
        AK16 = 10.0**CB = P(ALG)/P(ALL)
        '''
        CB = 5.70 - 15862 /sim.T
        self.EALL = 10.0**(-CB)

        '''
        Al(liq) + 0.5 O2(g) = AlO(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CC2 = - 2.43 + 13067 /sim.T
        AK17 = 10.0**CC = P(ALOG)/(P(ALL)*P(O2G)**0.5)
        HENCE P(ALOG) = 10.0**(CC-CB) * P(ALG) * P(O2G)**0.5
        '''
        CC = 3.04 - 2143 /sim.T
        self.EALOG = 10.0**(CC-CB)

        '''
        Al(liq) + O2(g) = AlO2(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CD2 = - 5.70 + 21159 /sim.T
        AK18 = 10.0**CD = P(ALO2G)/(P(ALL)*P(O2G))
        HENCE P(ALO2G) = 10.0**(CD-CB) * P(ALG) * P(O2G)
        '''
        CD = - 0.09 + 5.523 /sim.T
        self.EALO2G = 10.0**(CD-CB)

        '''
        2Al(liq) + 0.5 O2(g) = Al2O(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CE2 = - 9.32 + 41897 /sim.T
        AK19 = 10.0**CE = P(AL2OG)/P(ALL)**2 * P(O2G)**0.5
        HENCE P(AL2OG) = 10.0**(CE-2.0D0*CB) * P(ALG)**2 * P(O2G)**0.5
        '''
        CE = 2.04 + 10232 /sim.T
        self.EAL2OG = 10.0**(CE - 2.0 * CB)

        '''
        2Al(liq) + O2(g) = Al2O2(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CF2 = - 12.86 + 54600 /sim.T
        AK20 = 10.0**CF = P(AL2O2G)/(P(ALL)**2 * P(O2G))
        HENCE P(AL2O2G) = 10.0**(CF-2.0D0*CB) * P(ALG)**2 * PO2G
        '''
        CF = - 1.53 + 23021 /sim.T
        self.EAL2O2G = 10.0**(CF - 2.0 * CB)

        '''
        ######sim.Titanium ######
       sim.Ti(liq) + 0.5 O2(g) =sim.TiO(g)
        FOR P(TIG) INSTEAD OF P(TIL)    DC2 = - 3.33 + 23747 /sim.T
        AK22 = 10.0**DC = P(TIOG)/(P(TIL)*P(O2G)**0.5)
        HENCE P(TIL) = 10.0**(-DC) * P(TIOG) / P(O2G)**0.5
        '''
        DC = 4.31 - 2101 /sim.T
        self.ETIOG = 10.0**DC
        
        '''
       sim.Ti(liq) =sim.Ti(g)
        AK21 = 10.0**DB = P(TIG)/P(TIL)
        '''
        DB = 6.46 - 23025 /sim.T
        self.ETIL = 10.0**(-DB)
        
        '''
       sim.TiO2(liq) =sim.Ti(g) + 2O(g)
        JANAF 2nd ed. & supplements 1500-4000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99998
        AK20 = 10.0**DA = P(TIG)*P(OG)**2/P(TIO2L)
        HENCE P(TIO2L) = 10.0**(-DA) * P(TIG) * P(OG)**2
        '''
        DA = 21.07 - 95362 /sim.T
        self.ETIO2L = 10.0**(-DA)
    
        '''
        Ti(liq) + O2(g) = TiO2(g)
        FOR P(TIG) INSTEAD OF P(TIL)    DD2 = - 7.44 + 43028 /sim.T
        AK23 = 10.0**DD = P(TIO2G)/(P(TIL)*P(O2G))
        HENCE P(TIO2G) = 10.0**DD * P(TIL) * P(O2G)
        '''
        DD = - 0.41 + 17926 /sim.T
        self.ETIO2G = 10.0**DD

        '''
        ###### Sodium ###### 
        2Na(g) + O(g) = Na2O(liq)
        JANAF 2nd ed. & supplements
        correlation coefficient for linear fit = 
        '''
        EA = - 15.56 + 40286/sim.T
        self.ENA2OL = 10.0**EA

        '''
        Na(g) + O(g) = NaO(g)
        '''
        EB = - 1.43 + 1287/sim.T
        self.ENAOG = 10.0**(EB-E)

        '''
        2Na(g) = Na2(g)
        '''
        EC = - 4.31 + 4281/sim.T
        self.ENA2G = 10.0**EC

        '''
        2Na(g) + O(g) = Na2O(g)
        '''
        ED = - 7.00 + 11898/sim.T
        self.ENA2OG = 10.0**(ED-E)

        '''
        Na(g) = Na+(g) + e-(g)
        '''
        EE = 2.80 - 27851/sim.T
        self.ENACAT = 10**EE
        
        '''
        ####### Potassium ######
        2K(g) + O(g) = K2O(liq)
        '''
        FA = - 15.33 + 36735/sim.T
        self.EK2OL = 10.0**FA

        '''
        K(g) + O(g) = KO(g)
        '''
        FB = - 1.28 + 959/sim.T
        self.EKOG = 10.0**(FB-E)

        '''        
        2K(g) = K2(g)
        '''
        FC = - 3.94 + 2852/sim.T
        self.EK2G = 10.0**FC

        '''
        2K(g) + O(g) = K2O(g)
        '''
        FD = - 10.734 + 30817/sim.T
        self.EK2OG = 10.0**(FD-E)

        '''
        K(g) = K+(g) + e-(g)
        '''
        FE = 2.76 - 23760/sim.T
        self.EKCAT = 10**FE  

        '''
        ###### Thorium ######
        Th(liq) + O2(g) = ThO2(liq)
        '''
        GA = - 9.55 + 63948/sim.T
        self.ETHO2L = 10.0**GA

        '''
        Th(g) = Th(liq)
        '''
        GB = - 5.96 + 29600/sim.T
        self.ETHL = 10.0**GB
        
        '''
        Th(liq) + 0.5 O2(g) = ThO(g)
        '''
        GC = 2.75 + 3497/sim.T
        self.ETHOG = 10.0**GC
        
        '''
        Th(liq) + O2(g) = ThO2(g)
        '''
        GD = - 1.58 + 28875/sim.T
        self.ETHO2G = 10.0**GD

        '''
        ###### Uranium ######
        U(liq) + O2(g) = UO2(liq)
        '''
        HA = - 26.91 + 204359/sim.T
        self.EUO2L = 10.0**HA
        
        '''
        U(g) = U(liq)
        '''
        HB = - 5.75 + 25470/sim.T
        self.EUL = 10.0**HB
       
        '''
        U(liq) + 0.5 O2(g) = UO(g)
        '''
        HC = 3.02 + 1705/sim.T
        self.EUOG = 10.0**HC

        '''
        U(liq) + O2(g) = UO2(g)
        '''
        HD = - 1.19 + 26554/sim.T
        self.EUO2G = 10.0**HD

        '''
        U(liq) + 1.5 O2(g) = UO3(g)
        '''
        HE = - 4.24 +43710/sim.T
        self.EUO3G = 10.0**HE

        '''
        ###### Plutonium #####
        Pu(liq) + O2(g) = PuO2(liq)
        '''
        QA = - 29.86 + 200903/sim.T
        self.EPUO2L = 10.0**QA
       
        '''
        Pu(g) = Pu(liq)
        '''
        QB = - 4.79 + 17316/sim.T
        self.EPUL = 10.0**QB
       
        '''
        Pu(liq) + 0.5 O2(g) = PuO(g)
        '''
        QC = 2.40 + 6875/sim.T
        self.EPUOG = 10.0**QC
        
        '''
        Pu(liq) + O2(g) = PuO2(g)
        '''
        QD = - 1.76 + 25984/sim.T
        self.EPUO2G = 10.0**QD

    def melt_pressure_calculation(self,sim):
        '''
        Here the partial pressures of the oxides in the melt are calculated
        using the partial pressure of of the component gases and the melt
        vapor reaction. 
        '''
        ''' Si '''
        sim.presLiq['SiO2'] = self.ESIO2L * sim.presGas['SiO'] * sim.presGas['O2']**0.5
        sim.presLiq['Si']   = self.ESIL * sim.presGas['SiO'] * sim.presGas['O2']**(-0.5) 
        sim.presGas['Si']   = self.ESIG * sim.presLiq['Si'] 
        sim.presGas['O']    = self.EOG * sim.presGas['O2']**0.5
        # print(sim.presGas['O'])
        sim.presGas['SiO2'] = self.ESIO2G * sim.presLiq['Si'] * sim.presGas['O2']
      
        ''' Mg '''
        sim.presGas['Mg'] = self.EMGG * sim.presGas['MgO'] * sim.presGas['O2']**(-0.5)
        sim.presLiq['MgO'] = self.EMGOL *sim.presGas['Mg'] * sim.presGas['O']

        ''' Fe '''
        sim.presLiq['FeO'] = self.EFEOL * sim.presGas['Fe'] * sim.presGas['O']
        sim.presLiq['Fe']  = self.EFEL * sim.presGas['Fe']
        sim.presGas['FeO'] = self.EFEOG * sim.presGas['Fe'] * sim.presGas['O2']**0.5
        sim.presLiq['Fe2O3'] = self.EFE2O3L * sim.presGas['Fe']**2 * sim.presGas['O2']**1.5
        sim.presLiq['Fe3O4'] = self.EFE3O4L * sim.presGas['Fe']**3 * sim.presGas['O2']**2

        ''' Ca '''
        sim.presLiq['CaO'] = self.ECAOL * sim.presGas['Ca'] * sim.presGas['O']
        sim.presGas['CaO'] = self.ECAOG * sim.presGas['Ca'] * sim.presGas['O2']**0.5

        ''' Al '''
        sim.presLiq['Al2O3'] = self.EAL2O3L * sim.presGas['Al']**2 * sim.presGas['O']**3
        sim.presLiq['Al']    = self.EALL * sim.presGas['Al'] 
        sim.presGas['AlO']   = self.EALOG * sim.presGas['Al'] * sim.presGas['O2']**0.5
        sim.presGas['AlO2']  = self.EALO2G * sim.presGas['Al'] * sim.presGas['O2']
        sim.presGas['Al2O']  = self.EAL2OG * sim.presGas['Al']**2 * sim.presGas['O2']**0.5
        sim.presGas['Al2O2'] = self.EAL2O2G * sim.presGas['Al']**2 * sim.presGas['O2']

        ''' Ti ''' 
        sim.presLiq['Ti']   = self.ETIL * sim.presGas['Ti']
        sim.presGas['TiO']  = self.ETIOG * sim.presLiq['Ti'] * sim.presGas['O2']**0.5
        sim.presLiq['TiO2'] = self.ETIO2L * sim.presGas['Ti'] * sim.presGas['O']**2
        sim.presGas['TiO2'] = self.ETIO2G * sim.presLiq['Ti'] * sim.presGas['O2'] 

        ''' Na '''
        sim.presLiq['Na2O'] = self.ENA2OL * sim.presGas['Na']**2 * sim.presGas['O']
        sim.presGas['NaO']  = self.ENAOG * sim.presGas['Na'] * sim.presGas['O']
        sim.presGas['Na2']  = self.ENA2G * sim.presGas['Na']**2
        sim.presGas['Na2O'] = self.ENA2OG * sim.presGas['Na']**2 * sim.presGas['O']

        ''' K '''
        sim.presLiq['K2O'] = self.EK2OL * sim.presGas['K']**2 * sim.presGas['O']
        sim.presGas['KO'] = self.EKOG * sim.presGas['K'] * sim.presGas['O']
        sim.presGas['K2'] = self.EK2G * sim.presGas['K']**2
        sim.presGas['K2O'] = self.EK2OG * sim.presGas['K']**2 * sim.presGas['O']

        '''
        Ion chemistry:
        PENEG = P(e-,g), sim.presGas['NaCat'] = P(Na+,g), sim.presGas['KCat'] = P(K+,g)
        PENEG = sim.presGas['NaCat'] + sim.presGas['KCat']; 
        HENCE PENEG = DSQRT((sim.presGas['NaCat'] + PKCAT)*PENEG)
        '''
        # Energetic electron gas (???)
        sim.presGas['EnE'] = np.sqrt(self.ENACAT*sim.presGas['Na'] + self.EKCAT*sim.presGas['K'])

        if sim.presGas['Na'] != 0:
            sim.presGas['NaCat'] = self.ENACAT * sim.presGas['Na']/sim.presGas['EnE']
        else:
            sim.presGas['NaCat'] = 0

        if sim.presGas['K'] != 0:
            sim.presGas['KCat'] = self.EKCAT * sim.presGas['K']/sim.presGas['EnE']
        else: 
            sim.presGas['KCat'] = 0

    def number_density(self,sim):

        ''' Calculating for each species ''' 
        self.n = {}
        for gas in sim.presGas:
            self.n[gas] = self._pConv * sim.presGas[gas] / sim.T

        # # TODO find suitable names for these species
        # self.CENEG  = self._pConv * sim.presGas['EnE']  / sim.T
        # self.CNACAT = self._pConv * sim.presGas['NaCat'] / sim.T
        # self.CKCAT  = self._pConv * sim.presGas['KCat']  / sim.T

        ''' Calculating for each element ''' 
        self.n_el = {}
        self.n_el['Si'] = self.n['Si'] + self.n['SiO'] + self.n['SiO2']
        self.n_el['O']  = 2 * (self.n['O2'] + self.n['SiO2'] + self.n['AlO2'] + \
                                 self.n['Al2O2'] + self.n['TiO2']) + \
                            self.n['SiO'] + self.n['O'] + self.n['MgO'] + \
                            self.n['FeO'] + self.n['CaO'] + self.n['AlO'] + \
                            self.n['Al2O'] + self.n['TiO'] + self.n['NaO'] + \
                            self.n['Na2O'] + self.n['KO'] + self.n['K2O'] 
        self.n_el['Mg'] = self.n['MgO'] + self.n['Mg']
        self.n_el['Fe'] = self.n['FeO'] + self.n['Fe']
        self.n_el['Ca'] = self.n['CaO'] + self.n['Ca']
        self.n_el['Al'] = self.n['AlO'] + self.n['Al'] + self.n['AlO2'] + \
                            2 * (self.n['Al2O'] + self.n['Al2O2'])
        self.n_el['Ti'] = self.n['TiO'] + self.n['Ti'] + self.n['TiO2']
        self.n_el['Na'] = self.n['NaO'] + self.n['Na'] + self.n['NaCat'] + \
                            2 * (self.n['Na2'] + self.n['Na2O'])
        self.n_el['K']  = self.n['KO'] + self.n['K'] + self.n['KCat'] + \
                            2 * (self.n['K2O'] + self.n['K2'])

        ''' 
        Calculating the ratio of the number density of O in oxide gases to tot O
        number density
        '''
        self.oxideO_ratio = (self.n_el['Mg'] + self.n_el['Fe'] + self.n_el['Ca'] + \
                            2 * (self.n_el['Si'] + self.n_el['Ti']) + \
                            1.5 * (self.n_el['Al']) +  \
                            0.5 * (self.n_el['Na'] + self.n_el['K'])) / self.n_el['O']
        # print(self.n_el)

    def recompute_adjFact(self,sim):
        '''
        RECOMPUTE ADJUSTMENT FACTORS FOR THE KEY PRESSURES
        The activity of the oxides calculated during the activity
        calculations and from the gas chemistry must agree.
        
        If the oxide mole fraction for the element is zero, then there
        should not be any of that element in the vapor (A(element) = 0).
        
        The O2 abundance is governed by the most abundant oxide in the melt, 
        normally SiO2.  Once SiO2 is completely vaporized, AO2G is computed 
        from the remaining species in the melt, in order of volatility.
        '''

        # SiO
        if all(val != 0 for val in [sim.presLiq['SiO2'],sim.act_ox['SiO2']]):
            self.adjFact['SiO'] = 1 / (self.oxideO_ratio * (sim.presLiq['SiO2']\
                             / sim.act_ox['SiO2'])**0.5)
        elif sim.act_ox['SiO2'] == 0:
            self.adjFact['SiO'] = 0
        else:
            self.adjFact['SiO'] = 1

        # MgO
        if sim.presLiq['MgO'] != 0:
            self.adjFact['MgO'] = sim.act_ox['MgO'] / sim.presLiq['MgO'] 

        # Fe
        if sim.presLiq['FeO'] != 0 or sim.presLiq['Fe2O3'] != 0:
            self.adjFact['Fe'] = (sim.act_ox['FeO'] + sim.act_ox['Fe2O3'])\
                            / (sim.presLiq['FeO'] + sim.presLiq['Fe2O3'])
        elif sim.act_ox['FeO'] == 0:
            self.adjFact['Fe'] = 0
        else:
            self.adjFact['Fe'] = 1

        # Ca 
        if all(val != 0 for val in [sim.presLiq['CaO'],sim.act_ox['CaO']]):
            self.adjFact['Ca'] = (sim.act_ox['CaO'] / sim.presLiq['CaO'])**0.5
        elif sim.act_ox['CaO'] == 0:
            self.adjFact['Ca'] = 0
        else:
            self.adjFact['Ca'] = 1

        # Al
        if all(val != 0 for val in [sim.presLiq['Al2O3'],sim.act_ox['Al2O3']]):
            self.adjFact['Al'] = (sim.act_ox['Al2O3'] / sim.presLiq['Al2O3'])**0.5
        elif sim.act_ox['Al2O3'] == 0:
            self.adjFact['Al'] = 0
        else:
            self.adjFact['Al'] = 1

        # Ti
        if all(val != 0 for val in [sim.presLiq['TiO2'],sim.act_ox['TiO2']]):
            self.adjFact['Ti'] = (sim.act_ox['TiO2'] / sim.presLiq['TiO2'])**0.5
        elif sim.act_ox['TiO2'] == 0:
            self.adjFact['Ti'] = 0
        else:
            self.adjFact['Ti'] = 1

        # Na
        if all(val != 0 for val in [sim.presLiq['Na2O'],sim.act_ox['Na2O']]):
            self.adjFact['Na'] = (sim.act_ox['Na2O'] / sim.presLiq['Na2O'])**0.5
        elif sim.act_ox['Na2O'] == 0:
            self.adjFact['Na'] = 0
        else:
            self.adjFact['Na'] = 1

        # K
        if all(val != 0 for val in [sim.presLiq['K2O'],sim.act_ox['K2O']]):
            self.adjFact['K'] = (sim.act_ox['K2O']/sim.presLiq['K2O'])**0.5
        elif sim.act_ox['K2O'] == 0:
            self.adjFact['K'] = 0
        else:
            self.adjFact['K'] = 1

        '''
        Adjustment factor for oxygen is governed by the most abundant volatile
        metal oxide present in the melt.
        TODO: THIS ASSUMES THAT THE ABUNDANCES ARE ALWAYS ORDEREDAS IN THIS IF STATEMENT
        '''
        if sim.presLiq['SiO2'] != 0 and sim.act_ox['SiO2'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['SiO2']\
                            / sim.presLiq['SiO2']

        elif sim.act_ox['MgO'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['MgO']\
                            / sim.presLiq['MgO']

        elif sim.act_ox['FeO'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * (sim.act_ox['FeO']\
                            + sim.act_ox['Fe2O3'])\
                            / (sim.presLiq['FeO'] + sim.presLiq['Fe2O3'])

        elif sim.act_ox['CaO'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['CaO']\
                            / sim.presLiq['CaO']

        elif sim.act_ox['Al2O3'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['Al2O3']\
                            / sim.presLiq['Al2O3']

        elif sim.act_ox['TiO2'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['TiO2']\
                            / sim.presLiq['TiO2']

        elif sim.act_ox['Na2O'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['Na2O']\
                            / sim.presLiq['Na2O']

        elif sim.act_ox['K2O'] != 0:
            self.adjFact['O2'] = self.oxideO_ratio * sim.act_ox['K2O']\
                            / sim.presLiq['K2O']
        else:
            self.adjFact['O2'] = 1             


    # end number_density()

    def vapor_pressure_calculation(self,sim):

        iit = 0
        self.dif_range = 2.30359e-6 # Max amount that the adjFact can difer from 1

        while not all( 1-self.dif_range < fact < 1+self.dif_range or fact == 0 \
                   for fact in list(self.adjFact.values())) or iit == 0:

            # Adjust gas pressure according to adjustment factor
            sim.presGas = {gas : sim.presGas[gas]*self.adjFact[gas]\
                            for gas in sim._gasNames}

            # Calculate melt pressure using gas pressures
            self.melt_pressure_calculation(sim)

            # Calculate the number densities (TODO: this could go into vapor removal (?))
            self.number_density(sim)

            # Recompute the adjustment factors using the new partial pressures
            self.recompute_adjFact(sim)
            
            ''' While loop break '''
            iit += 1 # updating counter
            if iit >= 1e8: 
                raise RuntimeError('Max recursion limit reached while calculating adjustment factors.')

        # end while loop
    
    # end vapor_pressure_calculation()
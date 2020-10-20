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
        self.A = 22.13 - 94311.0 / self.T
         
        '''
        0.5O2 (g) = O(g)
        JANAF 2nd ed. 
        correlation coefficient for linear fit = 
        AK2 = 10.0**E = P(OG)/P(02G)**0.5

        HENCE 10.0**A = P(SIG)*(P(O2G)*10.0**2E)/P(SIO2L)
        HENCE P(SIO2L) = 10.0**(2.0*E-A)*P(SIG)*P(O2G)
        HENCE P(O2G) = 10.0**(-2.0*E) * P(OG)**2
        '''
        self.E = 3.47 - 13282 / self.T
        self.EOG = 10.0**self.E
        
        '''    
        Si(liq) = Si(g)
        AK3 = 10.0**B = P(SIG)/P(SIL)
        HENCE P(SIG) = 10.0**B*P(SIL)
        '''
        self.B = 6.00 - 20919 / self.T
        self.ESIG = 10.0**self.B

        '''
        Si(liq) + 0.5 O2(g) = SiO(g)
        FOR P(SIG) INSTEAD OF P(SIL)    C2 = - 3.67 + 29760 / T
        AK4 = 10.0**C = P(SIOG)/(P(SIL)*P(O2G)**0.5)
        HENCE P(SIL) = 10.0**-C*P(SIOG)/P(O2G)**0.5
        HENCE P(SIO2L) = 10.0**(2.0*E+B-A)*10.0**(-C)*P(SIOG)*P(O2G)**0.5
        HENCE P(SIO2L) = 10.0**(2.0*E+B-A-C)*P(SIOG)*P(O2G)**0.5
        '''
        self.C = 2.51 + 8207 / self.T
        self.ESIL = 10.0**(-self.C)
        self.ESIO2L = 10.0**(2.0*self.E + self.B - self.A - self.C)

        '''
        Si(liq) + O2(g) = SiO2(g)
        FOR P(SIG) INSTEAD OF P(SIL)    D2 = - 7.57 + 39595 / T
        AK5 = 10.0**D = P(SIO2G)/(P(SIL)*P(O2G))
        HENCE P(SIO2G) = 10.0**D*P(SIL)*P(O2G)
        HENCE P(SIO2G) = 10.0**D*10.0**-B*P(SIG)*P(O2G)
        HENCE P(SIO2G) = 10.0**(D-B) * P(SIG) * P(O2G)
        '''
        self.D = -1.44 + 18326 / self.T
        self.ESIO2G = 10.0**self.D

        '''
        ###### Magnesium ######
        MgO(liq) = Mg(g) + O(g)
        JANAF 2nd ed. & supplements 1500-5000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99994
        AK6 = 10.0**F = P(MGG)*P(OG)/P(MGOL)
        HENCE P(MGOL) = 10.0**(-F) * P(MGG) * P(OG)
        '''
        self.F = 12.56 - 46992 / self.T
        self.EMGOL = 10.0**(-self.F)
    
        '''
        Mg(g) + 0.5 O2(g) = MgO(g)
        AK8 = 10.0**G = P(MGOG)/(P(MGG)*P(O2G)**0.5)
        HENCE P(MGOG) = 10.0**G * P(MGG) * P(O2G)**0.5
        '''
        self.G = -1.19 + 3794 / self.T
        self.EMGG = 10.0**(-self.G)

        '''
        ###### Iron ######
        FeO(liq) = Fe(g) + O(g)
        JANAF 2nd ed. & supplements 2000-5000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99995
        AK9 = 10.0**AA = P(FEG)*P(OG)/P(FEOL)
        HENCE P(FEOL) = 10.0**(-AA) * P(FEG) * P(OG)
        '''
        self.AA = 12.06 - 44992 / self.T
        self.EFEOL = 10.0**(-self.AA)
        
        '''
        Fe(liq) = Fe(g)
        AK10 = 10.0**AB = P(FEG)/P(FEL)
        Fe(liq) + 0.5 O2(g) = FeO(g)
        FOR P(FEG) INSTEAD OF P(FEL)    AC2 = - 2.93 + 9945 / T
        AK11 = 10.0**AC = P(FEOG)/(P(FEL)*P(O2G)**0.5)
        HENCE P(FEOG) = 10.0**(AC-AB) * P(FEG) * P(O2G)**0.5
        '''
        self.AB = 6.35 -19704 / self.T
        self.EFEL = 10.0**(-self.AB)
        self.AC = 3.39 - 9951 / self.T
        self.EFEOG = 10.0**(self.AC-self.AB)

        '''
        2Fe (g) + 1.5 O2 (g) = Fe2O3 (liq)
        Fe2O3 (liq) cp data from IVTANTHERMO database (estimated)
        hematite enthalpy of fusion calculated from Sugawara & Akaogi 2004
        K = PFE2O3L / (PFEG**2 * PO2G**0.5)
        HENCE P(FE2O3L) = 10**AD * PFEG**2 * PO2G**0.5
        '''
        self.AD = -2.26722053113 + 7.56430936141329 / self.T
        self.EFE2O3L = 10**self.AD
        
        '''
        3Fe (g) + 2 O2 (g) = Fe3O4 (liq)
        Fe3O4 (liq) cp data from Barin 95
        magnetite enthalpy of fusion from JANAF 4th ed.
        K = PFE3O4L / (PFEG**3 * PO2G**0.5)
        HENCE P(FE3O4L) = 10**AE * PFEG**3 * PO2G**0.5
        '''
        self.AE = -3.19907301154 + 1.110526206139634 / self.T
        self.EFE3O4L = 10**self.AE

        '''
        ###### Calcium ######
        CaO(liq) = Ca(g) + O(g)
        JANAF 2nd ed. & supplements 2000-4500 K every 500 degrees 
        correlation coefficient for linear fit = -0.99998
        AK12 = 10.0**BA = P(CAG)*P(OG)/P(CAOL)
        HENCE P(CAOL) = 10.0**(-BA) * P(CAG) * P(OG)
        '''
        self.BA = 11.88 - 49586 / self.T
        self.ECAOL = 10.0**(-self.BA)
        
        '''
        Ca(g) + 0.5 O2(g) = CaO(g)
        AK14 = 10.0**BC = P(CAOG)/(P(CAG)*P(O2G)**0.5)
        HENCE P(CAOG) = 10.0**BC * P(CAG) * P(O2G)**0.5
        '''
        self.BC = -1.61 + 6128 / self.T
        self.ECAOG = 10.0**self.BC

        '''
        ###### Aluminum ###### 

        Al2O3(liq) = 2Al(g) + 3O(g)
        JANAF supplements 1500-4000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99995
        AK15 = 10.0**CA = P(ALG)**2*P(OG)**3/P(AL2O3L)
        HENCE P(AL2O3L) = 10.0**(-CA) * P(ALG)**2 * P(OG)**3
        '''
        self.CA = 35.83 - 153255 / self.T
        self.EAL2O3L = 10.0**(-self.CA)
    
        '''
        Al(liq) = Al(g)
        AK16 = 10.0**CB = P(ALG)/P(ALL)
        '''
        self.CB = 5.70 - 15862 / self.T
        self.EALL = 10.0**(-self.CB)

        '''
        Al(liq) + 0.5 O2(g) = AlO(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CC2 = - 2.43 + 13067 / T
        AK17 = 10.0**CC = P(ALOG)/(P(ALL)*P(O2G)**0.5)
        HENCE P(ALOG) = 10.0**(CC-CB) * P(ALG) * P(O2G)**0.5
        '''
        self.CC = 3.04 - 2143 / self.T
        self.EALOG = 10.0**(self.CC-self.CB)

        '''
        Al(liq) + O2(g) = AlO2(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CD2 = - 5.70 + 21159 / T
        AK18 = 10.0**CD = P(ALO2G)/(P(ALL)*P(O2G))
        HENCE P(ALO2G) = 10.0**(CD-CB) * P(ALG) * P(O2G)
        '''
        self.CD = - 0.09 + 5.523 / self.T
        self.EALO2G = 10.0**(self.CD-self.CB)

        '''
        2Al(liq) + 0.5 O2(g) = Al2O(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CE2 = - 9.32 + 41897 / T
        AK19 = 10.0**CE = P(AL2OG)/P(ALL)**2 * P(O2G)**0.5
        HENCE P(AL2OG) = 10.0**(CE-2.0D0*CB) * P(ALG)**2 * P(O2G)**0.5
        '''
        self.CE = 2.04 + 10232 / self.T
        self.EAL2OG = 10.0**(self.CE - 2.0 * self.CB)

        '''
        2Al(liq) + O2(g) = Al2O2(g)
        FOR P(ALG) INSTEAD OF P(ALL)    CF2 = - 12.86 + 54600 / T
        AK20 = 10.0**CF = P(AL2O2G)/(P(ALL)**2 * P(O2G))
        HENCE P(AL2O2G) = 10.0**(CF-2.0D0*CB) * P(ALG)**2 * PO2G
        '''
        self.CF = - 1.53 + 23021 / self.T
        self.EAL2O2G = 10.0**(self.CF - 2.0 * self.CB)

        '''
        ###### Titanium ######
        Ti(liq) + 0.5 O2(g) = TiO(g)
        FOR P(TIG) INSTEAD OF P(TIL)    DC2 = - 3.33 + 23747 / T
        AK22 = 10.0**DC = P(TIOG)/(P(TIL)*P(O2G)**0.5)
        HENCE P(TIL) = 10.0**(-DC) * P(TIOG) / P(O2G)**0.5
        '''
        self.DC = 4.31 - 2101 / self.T
        self.ETIOG = 10.0**self.DC
        
        '''
        Ti(liq) = Ti(g)
        AK21 = 10.0**DB = P(TIG)/P(TIL)
        '''
        self.DB = 6.46 - 23025 / self.T
        self.ETIL = 10.0**(-self.DB)
        
        '''
        TiO2(liq) = Ti(g) + 2O(g)
        JANAF 2nd ed. & supplements 1500-4000 K every 500 degrees 
        correlation coefficient for linear fit = -0.99998
        AK20 = 10.0**DA = P(TIG)*P(OG)**2/P(TIO2L)
        HENCE P(TIO2L) = 10.0**(-DA) * P(TIG) * P(OG)**2
        '''
        self.DA = 21.07 - 95362 / self.T
        self.ETIO2L = 10.0**(-self.DA)
    
        '''
        Ti(liq) + O2(g) = TiO2(g)
        FOR P(TIG) INSTEAD OF P(TIL)    DD2 = - 7.44 + 43028 / T
        AK23 = 10.0**DD = P(TIO2G)/(P(TIL)*P(O2G))
        HENCE P(TIO2G) = 10.0**DD * P(TIL) * P(O2G)
        '''
        self.DD = - 0.41 + 17926 / self.T
        self.ETIO2G = 10.0**self.DD

        '''
        ###### Sodium ###### 
        2Na(g) + O(g) = Na2O(liq)
        JANAF 2nd ed. & supplements
        correlation coefficient for linear fit = 
        '''
        self.EA = - 15.56 + 40286/self.T
        self.ENA2OL = 10.0**self.EA

        '''
        Na(g) + O(g) = NaO(g)
        '''
        self.EB = - 1.43 + 1287/self.T
        self.ENAOG = 10.0**(self.EB-self.E)

        '''
        2Na(g) = Na2(g)
        '''
        self.EC = - 4.31 + 4281/self.T
        self.ENA2G = 10.0**self.EC

        '''
        2Na(g) + O(g) = Na2O(g)
        '''
        self.ED = - 7.00 + 11898/self.T
        self.ENA2OG = 10.0**(self.ED-self.E)

        '''
        Na(g) = Na+(g) + e-(g)
        '''
        self.EE = 2.80 - 27851/self.T
        self.ENACAT = 10**self.EE
        
        '''
        ####### Potassium ######
        2K(g) + O(g) = K2O(liq)
        '''
        self.FA = - 15.33 + 36735/self.T
        self.EK2OL = 10.0**self.FA

        '''
        K(g) + O(g) = KO(g)
        '''
        self.FB = - 1.28 + 959/self.T
        self.EKOG = 10.0**(self.FB-self.E)

        '''        
        2K(g) = K2(g)
        '''
        self.FC = - 3.94 + 2852/self.T
        self.EK2G = 10.0**self.FC

        '''
        2K(g) + O(g) = K2O(g)
        '''
        self.FD = - 10.734 + 30817/self.T
        self.EK2OG = 10.0**(self.FD-self.E)

        '''
        K(g) = K+(g) + e-(g)
        '''
        self.FE = 2.76 - 23760/self.T
        self.EKCAT = 10**self.FE  

        '''
        ###### Thorium ######
        Th(liq) + O2(g) = ThO2(liq)
        '''
        self.GA = - 9.55 + 63948/self.T
        self.ETHO2L = 10.0**self.GA

        '''
        Th(g) = Th(liq)
        '''
        self.GB = - 5.96 + 29600/self.T
        self.ETHL = 10.0**self.GB
        
        '''
        Th(liq) + 0.5 O2(g) = ThO(g)
        '''
        self.GC = 2.75 + 3497/self.T
        self.ETHOG = 10.0**self.GC
        
        '''
        Th(liq) + O2(g) = ThO2(g)
        '''
        self.GD = - 1.58 + 28875/self.T
        self.ETHO2G = 10.0**self.GD

        '''
        ###### Uranium ######
        U(liq) + O2(g) = UO2(liq)
        '''
        self.HA = - 26.91 + 204359/self.T
        self.EUO2L = 10.0**self.HA
        
        '''
        U(g) = U(liq)
        '''
        self.HB = - 5.75 + 25470/self.T
        self.EUL = 10.0**self.B
       
        '''
        U(liq) + 0.5 O2(g) = UO(g)
        '''
        self.HC = 3.02 + 1705/self.T
        self.EUOG = 10.0**self.HC

        '''
        U(liq) + O2(g) = UO2(g)
        '''
        self.HD = - 1.19 + 26554/self.T
        self.EUO2G = 10.0**self.HD

        '''
        U(liq) + 1.5 O2(g) = UO3(g)
        '''
        self.HE = - 4.24 +43710/self.T
        self.EUO3G = 10.0**self.HE

        '''
        ###### Plutonium #####
        Pu(liq) + O2(g) = PuO2(liq)
        '''
        self.QA = - 29.86 + 200903/self.T
        self.EPUO2L = 10.0**self.QA
       
        '''
        Pu(g) = Pu(liq)
        '''
        self.QB = - 4.79 + 17316/self.T
        self.EPUL = 10.0**self.QB
       
        '''
        Pu(liq) + 0.5 O2(g) = PuO(g)
        '''
        self.QC = 2.40 + 6875/self.T
        self.EPUOG = 10.0**self.QC
        
        '''
        Pu(liq) + O2(g) = PuO2(g)
        '''
        self.QD = - 1.76 + 25984/self.T
        self.EPUO2G = 10.0**self.QD


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

    def ion_chemistry(self, presGas, presLiq):

        '''
        Compute the partial pressures of the vapor species and activities of the oxides.
        
        Ion chemistry:
            PENEG = P(e-,g), PNACAT = P(Na+,g), PKCAT = P(K+,g)
            PENEG = PNACAT + PKCAT; 
            HENCE PENEG = DSQRT((PNACAT + PKCAT)*PENEG)
        '''

        # TODO: Look into placing all the relevant variables in thermo data
        ''' Si '''
        presLiq['SiO2'] = self.ESIO2L * presGas['SiO']**0.5
        presLiq['Si']   = self.ESIL * presGas['SiO'] * presGas['O2']**(-0.5) 
        presGas['Si']   = self.ESIG * presLiq['Si'] 
        presGas['O']    = self.EOG * presGas['O2']**0.5
        presGas['SiO2'] = self.ESIO2G * presLiq['Si'] * presGas['O2']

        # PSIO2L = ESIO2L * PSIOG * DSQRT(PO2G)
        # PSIL = ESIL * PSIOG * PO2G**(-0.5)
        # PSIG = ESIG * PSIL
        # POG = EOG * DSQRT(PO2G)
        # PSIO2G = ESIO2G * PSIL * PO2G

        ''' Mg '''
        presGas['Mg'] = self.EMGG * presGas['MgO'] * presGas['O2']**(-0.5)
        presLiq['Mg'] = self.EMGOL *presGas['Mg'] * presGas['O']

        # PMGG = EMGG * PMGOG * PO2G**(-0.5)
        # PMGOL = EMGOL * PMGG * POG

        ''' Fe '''
        presLiq['FeO'] = self.EFEOL * presGas['Fe'] * presGas['O']
        presLiq['Fe']  = self.EFEL * presGas['Fe']
        presGas['FeO'] = self.EFEOG * presGas['Fe'] * presGas['O2']**0.5
        presLiq['Fe2O3'] = self.EFE2O3L * presGas['Fe']**2 * presGas['O2']**1.5
        presLiq['Fe3O4'] = self.EFE3O4L * presGas['Fe']**3 * presGas['O2']**2

        # PFEOL = EFEOL * PFEG * POG
        # PFEL = EFEL * PFEG
        # PFEOG = EFEOG * PFEG * DSQRT(PO2G)
        # PFE2O3L = EFE2O3L * PFEG**2 * PO2G**1.5
        # PFE3O4L = EFE3O4L * PFEG**3 * PO2G**2

        ''' Ca '''
        presLiq['CaO'] = self.ECAOL * presGas['Ca'] * presGas['O']
        presGas['CaO'] = self.ECAOG * presGas['Ca'] * presGas['O2']

        # PCAOL = ECAOL * PCAG * POG
        # PCAOG = ECAOG * PCAG * DSQRT(PO2G)

        ''' Al '''
        presLiq['Al2O3'] = self.EAL2O3L * presGas['Al']**2 * presGas['O']**3
        presLiq['Al']    = self.EALL * presGas['Al'] 
        presGas['AlO']   = self.EALOG * presGas['Al'] * presGas['O2']**0.5
        presGas['AlO2']  = self.EALO2G * presGas['Al'] * presGas['O2']
        presGas['Al2O']  = self.EAL2OG * presGas['Al']**2 * presGas['O2']**0.5
        presGas['Al2O2'] = self.EAL2O2G * presGas['Al']**2 * presGas['O2']

        # PAL2O3L = EAL2O3L * PALG**2 * POG**3
        # PALL = EALL * PALG
        # PALOG = EALOG * PALG * DSQRT(PO2G)
        # PALO2G = EALO2G * PALG * PO2G
        # PAL2OG = EAL2OG * PALG**2 * DSQRT(PO2G)
        # PAL2O2G = EAL2O2G * PALG**2 * PO2G

        ''' Ti ''' 
        presLiq['Ti']   = self.ETIL * presGas['Ti']
        presGas['TiO']  = self.ETIOG * presLiq['Ti'] * presGas['O2']**0.5
        presLiq['TiO2'] = self.ETIO2L * presGas['Ti'] * presGas['O']**2
        presGas['TiO2'] = self.ETIO2G * presLiq['Ti'] * presGas['O2'] 

        # PTIL = ETIL * PTIG
        # PTIOG = ETIOG * PTIL * DSQRT(PO2G)
        # PTIO2L = ETIO2L * PTIG * POG**2
        # PTIO2G = ETIO2G * PTIL * PO2G

        ''' Na '''
        presLiq['Na2O'] = self.ENA2OL * presGas['Na']**2 * presGas['O']
        presGas['NaO']  = self.ENAOG * presGas['Na'] * presGas['O']
        presGas['Na2']  = self.ENA2G * presGas['Na']**2
        presGas['Na2O'] = self.ENA2OG * presGas['Na']**2 * presGas['O']

        # PNA2OL = ENA2OL * PNAG**2 * POG
        # PNAOG = ENAOG * PNAG * POG
        # PNA2G = ENA2G * PNAG**2
        # PNA2OG = ENA2OG * PNAG**2 * POG

        ''' K '''
        presLiq['K2O'] = self.EK2OL * presGas['K']**2 * presGas['O']
        presGas['KO'] = self.EKOG * presGas['K'] * presGas['O']
        presGas['K2'] = self.EK2G * presGas['K']**2
        presGas['K2O'] = self.EK2OG * presGas['K']**2 * presGas['O']

        # PK2OL = EK2OL * PKG**2 * POG
        # PKOG = EKOG * PKG * POG
        # PK2G = EK2G * PKG**2
        # PK2OG = EK2OG * PKG**2 * POG

        # TODO: findin out what ENACAT, EKCAT and PENEG stand for
        self.PENEG = np.sqrt(self.ENACAT*presGas['Na'] + self.EKCAT*presGas['K'])

        # PENEG = DSQRT(ENACAT*PNAG + EKCAT*PKG)

        if presGas['Na'] != 0:
            self.PNACAT = self.ENACAT * presGas['Na']/self.PENEG
        else:
            self.PNACAT = 0

        # IF(PNAG .NE. 0.0D0) THEN
        # PNACAT = ENACAT * PNAG / PENEG
        # ELSE
        # PNACAT = 0.0D0
        # ENDIF

        if presGas['K'] != 0:
            self.PKCAT = self.EKCAT * presGas['K']/self.PENEG
        else: 
            self.PKCAT = 0

        # IF(PKG .NE. 0.0D0) THEN
        # PKCAT = EKCAT * PKG / PENEG
        # ELSE
        # PKCAT = 0.0D0
        # ENDIF

    # end ion_chemistry()
        
# end melt_activities():
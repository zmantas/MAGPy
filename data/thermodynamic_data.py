'''
Gas chemistry thermodynamic data
'''

class melt_activities():
    '''
    Thermodynamic data for activities in the melt
    '''
    def __init__(self,T):
        '''
        2MgO(liq) + SiO2(liq) = Mg2SiO4(liq)
        log10K(Mg2SiO4) = - 34.08 + 141582/T
        -2log10K(MgO)   =  18.10 - 67242/T
         -log10K(SiO2)  =  15.04 - 66906/T
        '''
        self.AKMG2 = 10**(- 0.94 + 7434/T)

        '''
        MgO(liq) + SiO2(liq) = MgSiO3(liq)
        log10K(MgSiO3) = - 23.67 + 102856/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.AKMG1 = 10**( 0.42 + 2329/T)

        '''
        MgO(liq) + Al2O3(liq) = MgAl2O4(liq)
        log10K(MgAl2O4) = - 31.55 + 142219/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.AKMG3 = 10**(1.18 + 464/T)

        '''
        MgO(liq) + TiO2(liq) = MgTiO3(liq)
        log10K(MgTiO3) = - 22.54 + 103180/T
        -log10K(MgO) = 9.05 - 33621/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.AKMG4 = 10**(- 0.13 + 3246/T)

        '''
        MgO(liq) + 2TiO2(liq) = MgTi2O5(liq)
        log10K(MgTi2O5) = - 35.26 + 169092/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10(TiO2) = 26.72 - 132626/T
        '''
        self.AKMG5 = 10**(0.51 + 2845/T)

        '''
        2MgO(liq) + TiO2(liq) = Mg2TiO4(liq)
        log10K(Mg2TiO4) = - 30.79 + 137367/T
        -2log10K(MgO) = 18.10 - 67242/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.AKMG6 = 10**(0.67 + 3812/T)

        '''
        2MgO(liq) + 2Al2O3(liq) + 5SiO2(liq) = Mg2Al4Si5O18(liq)
        log10K(Mg2Al4Si5O18) = - 132.38 + 618040/T
        -2log10K(MgO) = 18.10 - 67242/T
        -2log10K(Al2O3) = 47.36 - 216268/T
        -5log10K(SiO2) = 75.20 - 334530/T
        '''
        self.AKMG7 = 10**7.48

        '''
        3Al2O3(liq) + 2SiO2(liq) = Al6Si2O13(liq)
        log10K(Al6Si2O13) = - 104.06 + 467589/T
        -3log10K(Al2O3) = 71.04 - 324402/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.AKAL1 = 10**(- 2.94 + 9375/T)

        '''
        CaO(liq) + Al2O3(liq) = CaAl2O4(liq)
        log10K(CaAl2O4) = - 33.93 + 154384/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.AKCA1 = 10**(- 1.89 + 10060/T)

        '''
        CaO(liq) + 2Al2O3(liq) = CaAl4O7(liq)
        log10K(CaAl4O7) = - 56.31 + 262171/T
        -log10K(CaO) = 8.36 - 36190/T
        -2log10K(Al2O3) = 47.36 - 216268/T
        '''
        self.AKCA2 = 10**(- 0.59 + 9713/T)

        '''
        12CaO(liq) + 7Al2O3(liq) = Ca12Al14O33(liq)
        log10K(Ca12Al14O33) = - 272.38 + 1263457/T
        -12log10K(CaO) = 100.32 - 434280/T
        -7log10K(Al2O3) = 165.76 - 756938/T
        '''
        self.AKCA3 = 10**(-6.30 + 72239/T)

        '''
        CaO(liq) + SiO2(liq) = CaSiO3(liq)
        log10K(CaSiO3) = - 22.86 + 108664/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.AKCA4 = 10**(0.54 + 5568/T)

        '''
        CaO(liq) + Al2O3(liq) + 2SiO2(liq) = CaAl2Si2O8(liq)
        log10K(CaAl2Si2O8) = - 59.49 + 283462/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(Al2O3) = 23.68 - 108134/T
        -2log10(SiO2) = 30.08 - 133812/T
        '''
        self.AKCA5 = 10**(+ 2.63 + 5326/T)

        '''
        CaO(liq) + MgO(liq) + 2SiO2(liq) = CaMgSi2O6(liq)
        log10K(CaMgSi2O6) = -46.03 + 212108/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.AKCA6 = 10**(1.46 + 8485/T)

        '''
        2CaO(liq) + MgO(liq) + 2SiO2(liq) = Ca2MgSi2O7(liq)
        log10K(Ca2MgSi2O7) = - 55.22 + 255140/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(MgO) = 9.05 - 33621/T
        -2log10K(SiO2) = 30.08 - 133812/T
        '''
        self.AKCA7 = 10**(0.63 + 15327/T)

        '''
        2CaO(liq) + Al2O3(liq) + SiO2(liq) = Ca2Al2SiO7(liq)
        log10K(Ca2Al2SiO7) = - 53.43 + 258130/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(Al2O3) = 23.68 - 108134/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.AKCA8 = 10**(2.01 + 10710/T)

        '''
        CaO(liq) + TiO2(liq) = CaTiO3(liq)
        log10K(CaTiO3) = - 21.80 + 109558/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.AKCA9 = 10**(- 0.08 + 7055/T)

        '''
        2CaO(liq) + SiO2(liq) = Ca2SiO4(liq)
        log10K(Ca2SiO4) = - 31.13 + 147702/T
        -2log10K(CaO) = 16.72 - 72380/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.AKCA10 = 10**(0.63 + 8416/T)

        '''
        CaO(liq) + TiO2(liq) + SiO2(liq) = CaTiSiO5(liq)
        log10K(CaTiSiO5) = - 36.94 + 179480/T
        -log10K(CaO) = 8.36 - 36190/T
        -log10K(TiO2) = 13.36 - 66313/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.AKCA11 = 10**(- 0.18 + 10071/T)

        '''
        CaO(liq) + 6Al2O3(liq) = CaAl12O19(liq)
        log10K(CaAl12O19) = -154.23 + 707606/T
        -log10K(CaO) = 8.36 - 36190/T
        -6log10K(Al2O3) = 142.08 - 648804/T
        '''
        self.AKCA12 = 10**(- 3.79 + 22612/T)

        '''
        FeO(liq) + TiO2(liq) = FeTiO3(liq)
        log10K(FeTiO3) = - 22.14 + 100392/T
        -log10K(FeO) = 8.27 - 30510/T
        -log10K(TiO2) = 13.36 - 66313/T
        '''
        self.AKFE1 = 10**(- 0.51 + 3569/T)

        '''
        2FeO(liq) + SiO2(liq) = Fe2SiO4(liq)
        log10K(Fe2SiO4) = - 32.21 + 131029/T
        -2log10K(FeO) = 16.54 - 61020/T
        -log10K(SiO2) = 15.04 - 66906/T
        '''
        self.AKFE2 = 10**(- 0.63 + 3103/T)

        '''
        FeO(liq) + Al2O3(liq) = FeAl2O4(liq)
        log10K(FeAl2O4) = - 33.71 + 144336/T
        -log10K(FeO) = 8.27 - 30510/T
        -log10K(Al2O3) = 23.68 - 108134/T
        '''
        self.AKFE3 = 10**(- 1.76 + 5692/T)

        '''
        FeO (liq) + Fe2O3 (liq) = Fe3O4 (liq)
        Fe3O4 data from Barin 1995
        '''
        self.AKFE4 = 10**(-4.385894544D-1 + 4.3038155175436D03 / T -
             * 3.1050205223386055D6 / T**2.0D0)

        '''
        Na2O(liq) + SiO2(liq) = Na2SiO3(liq)
        '''
        self.AKNA1 = 10**(- 1.33 + 13870/T)

        '''
        Na2O(liq) + 2SiO2(liq) = Na2Si2O5(liq)
        '''
        self.AKNA2 = 10**(- 1.39 + 15350/T)

        ''' 
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + SiO2(liq) = NaAlSiO4(liq)
        '''
        self.AKNA3 = 10**(0.65 + 6997/T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) + 3SiO2(liq) = NaAlSi3O8(liq)
        '''
        self.AKNA4 = 10**(1.29 + 8788/T)

        '''
        0.5 Na2O(liq) + 0.5 Al2O3(liq) = NaAlO2(liq)
        '''
        self.AKNA5 = 10**(0.55 + 3058/T)

        '''
        Na2O(liq) + TiO2(liq) = Na2TiO3(liq)
        '''
        self.AKNA6 = 10**(- 1.38 + 15445/T)

        '''
        0.5 NA2O(liq) + 0.5 Al2O3(liq) + 2SiO2(liq) = NAAlSi2O6(liq)
        '''
        self.AKNA7 = 10**(- 1.02 + 9607/T)

        '''
        K2O(liq) + SiO2(liq) = K2SiO3(liq)
        '''
        self.AKK1 = 10**(0.2692 + 12735/T)

        '''
        K2O(liq) + 2SiO2(liq) = K2Si2O5(liq)
        '''
        self.AKK2 = 10**(0.3462 + 14685/T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + SiO2(liq) = KAlSiO4(liq)
        '''
        self.AKK3 = 10**(0.97 + 8675/T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + 3SiO2(liq) = KAlSi3O8(liq)
        '''
        self.AKK4 = 10**(1.11 + 11229/T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) = KAlO2(liq)
        '''
        self.AKK5 = 10**(0.72 + 4679/T)

        '''
        0.5 K2O(liq) + 0.5 Al2O3(liq) + 2SiO2(liq) = KAlSi2O6(liq)
        '''
        self.AKK6 = 10**(1.53 + 10125/T)

        '''
        K2O(liq) + 4SiO2 (liq) = K2Si4O9 (liq)
        '''
        self.AKK7 = 10**(-0.9648 + 17572 / T)

        '''
        0.5K2O(liq) + CaO(liq) + 0.5Al2O3(liq) + 2SiO2(liq)=KCaAlSi2O7(liq)
        '''
        self.AKK8 = 10**(4.2983 + 17037/T)

    # __init__()
        
# end melt_activities():
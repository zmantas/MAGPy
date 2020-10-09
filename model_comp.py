import numpy as np

class model_comp():
    '''
    Imports and contains all relevant information about the initial composition.
    '''
    def __init__(self,input_fname):

        self.input_fname = input_fname

        ''' Defining oxide names '''
        self._oxide_names = "SiO2",\
                            "MgO",\
                            "Al2O3",\
                            "TiO2",\
                            "Fe2O3",\
                            "FeO",\
                            "CaO",\
                            "Na2O",\
                            "K2O"

        ''' Reading data from input file '''
        print(f'Reading input data from: {input_fname}')
        model = np.genfromtxt(input_fname,dtype=None,delimiter='\n',encoding='UTF-8',skip_header=2)
        
        # Assigning wt values to corresponding oxide names
        self.wt_oxides = {}
        for i,name in enumerate(self._oxide_names):
            self.wt_oxides[name] = np.float(model[i])

    # end __init__()

    def print(self):
 
        ''' Prints composition of the model '''
        print(f'Model read from {self.input_fname} is composed of:')
        for name, wt in self.wt_oxides.items():
            print(f'{name}: {wt}')

    # end print()

# end model_comp()
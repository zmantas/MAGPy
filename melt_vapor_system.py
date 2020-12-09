# Standard libraries
import numpy as np

class system():
	'''
	Class containing all relevant information for the melt_vapor system.
	'''

	def __init__(self,input_fname):

		'''Constants'''
		self._avog = 6.023e23 # Avogadro's number

		''' 
		Defining variables
		'''
		# Oxide names
		self._oxideNames = 'SiO2','MgO','Al2O3','TiO2','Fe2O3','FeO','CaO',\
						   'Na2O','K2O'
		# Metal names 
		self._metalNames = 'Si','Mg','Al','Ti','Fe3','Fe','Ca','Na','K'

		# Metal to oxide dictionary
		self._metal2oxide = {self._metalNames[i] : self._oxideNames[i] \
							 for i in range(len(self._metalNames))}

		# Oxide activities
		self.act_ox = {ox : 0 for ox in self._oxideNames}

		'''
		Import standard data values
		'''

		# File names
		mwOxides_fname = 'data/weights_oxides.csv'
		wMetals_fname = 'data/weights_metals.csv'

		# Importing molecular oxide and metal weights
		self._mwOxides = dict(np.genfromtxt(mwOxides_fname,delimiter= ',',\
							  names=True,skip_header=1,dtype=None,\
							  encoding = 'UTF-8'))
		self._wMetals = dict(np.genfromtxt(wMetals_fname,delimiter= ',',\
							 names=True,skip_header=1,dtype=None,\
							 encoding = 'UTF-8',autostrip=True))

		# Evaluating the multiplications and dividing by _avog 
		# (mol wt. in g/mole)
		self._wMetals.update((x, eval(y)/self._avog) \
							  for x, y in self._wMetals.items())
		

		'''
		Importing composition

		Note: Names must be given for each oxide.

		TODO: Build check input data and put in other function

		'''        
		self.comp_init = dict(np.genfromtxt(input_fname,encoding='UTF-8',
									   skip_header=2,dtype=None,delimiter=','))
		self.totWt = sum(self.comp_init.values())


		# Number of moles of each oxide
		self.molOx = {ox : self.comp_init[ox] / self._mwOxides[ox] if\
					  self.comp_init[ox] != 0 else 0 for ox in self._mwOxides}
		self.totMol = sum(self.molOx.values())

		# Mole percentage of each oxide
		self.molOx_perc = {ox : 100 * self.molOx[ox] / self.totMol if\
						   self.molOx != 0 else 0 for ox in self.molOx}
		self.totPerc = sum(self.molOx_perc.values()) 

		'''
		------ Calculating elemental abundances ------
		
		- If there is Si in the composition, the abundances are normalised to
		  Si = 1e6, else avog constant is used. 

		- Equation: 
			abundance = # of metal atoms in oxide * moles * 1e6/moles of SiO2
		'''
		self.abEl = {}
		if self.molOx['SiO2'] != 0:            
			self.abEl['Si'] = 1e6
			self.abEl['Mg'] = self.molOx['MgO']  * 1e6 / self.molOx['SiO2']
			self.abEl['Al'] = self.molOx['Al2O3'] * 2 * 1e6 / self.molOx['SiO2']
			self.abEl['Ti'] = self.molOx['TiO2'] * 1e6 / self.molOx['SiO2']
			self.abEl['Fe'] = (self.molOx['FeO'] + 2 * self.molOx['Fe2O3'])\
							  * 1e6 / self.molOx['SiO2']
			self.abEl['Ca'] = self.molOx['CaO']  * 1e6 / self.molOx['SiO2']
			self.abEl['Na'] = self.molOx['Na2O'] * 2 * 1e6 / self.molOx['SiO2']
			self.abEl['K']  = self.molOx['K2O']  * 2 * 1e6 / self.molOx['SiO2']

		else:
			self.abEl['Si'] = self.molOx['SiO2'] * self._avog
			self.abEl['Mg'] = self.molOx['MgO']  * self._avog
			self.abEl['Al'] = self.molOx['Al2O3'] * 2 * self._avog
			self.abEl['Ti'] = self.molOx['TiO2'] * self._avog
			self.abEl['Fe'] = (self.molOx['FeO'] + 2 * self.molOx['Fe2O3']) \
							  * self._avog
			self.abEl['Ca'] = self.molOx['CaO']  * self._avog
			self.abEl['Na'] = self.molOx['Na2O'] * 2 * self._avog
			self.abEl['K']  = self.molOx['K2O']  * 2 * self._avog 


		''' Renomarlizing the abundances ''' 
		# Total atomic abundance of all the elements (except 0)
		self.abETot = sum(self.abEl.values())
		# Molecular abundance of all the oxides
		self.abETot_ox = self.abEl['Si'] + self.abEl['Mg'] + self.abEl['Fe'] + \
				  self.abEl['Ca'] + self.abEl['Ti'] + \
				  0.5 * (self.abEl['Al'] + self.abEl['Na'] + self.abEl['K']) 

		# Relative abundance of the metals per oxide
		self.fAbOx = {el : 0.5 * self.abEl[el] / self.abETot_ox if \
							el in ['Al','Na','K'] else self.abEl[el] /\
							self.abETot_ox for el in self.abEl} 
		
		# Relative abundance of the metals per atom							
		self.fAbAtom = {el : self.abEl[el] / self.abETot_ox for el in self.abEl} 

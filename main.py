# Standard libraries
import sys

# Local modules
from melt_vapor_system import system
from melt_activity import melt_activity


def main():

	''' 
	Setting initial values 
	'''

	# Temperature for which the calculations will be done '''
	T = 2200 # Temperature of magma in Kelvin

	# File names
	# input_fname = 'input/BSE-initial.dat'
	input_fname = 'input/ic_Komatiite.dat'
	output_fname = 'output/MAGMA.OUT'

	'''
	Initialising classes
	'''
	sim = system(input_fname)
	melt = melt_activity(T,sim)


	'''
	Running calculations
	'''
	melt.melt_activity_calculation(sim)



if __name__ == "__main__":
	sys.exit(main())
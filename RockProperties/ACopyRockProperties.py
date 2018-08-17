# My Python rock properties modules

import numpy as np
import scipy as sp
import MyModule as mm
from math import sqrt

def get_Vp_from_K_Mu_Rho(K, Mu, Rho):
	"""
	Computes the compressional velocity from bulk modulus, shear modulus and density.
	The units must be as described below.
	
		Input Variables:
			K   = Bulk modulus           (Gpa)
			Mu  = Shear Modulus          (Gpa)
			Rho = Density                (g/cc)
	
		Output Variable:
			Vp  = Compressional Velocity (m/s)
	
	"""
	# Conversion to (mks) system for computation
	K = 1000000000*K         # Conversion from GPa to Pa
	Mu = 1000000000*Mu       # Conversion from GPa to Pa
	Rho = 1000 * Rho         # Conversion from g/(cm^3) to Kg/(m^3)
	
	Vp = sqrt( (K + 4*Mu/3)/Rho )
	
	
	return Vp
	
def get_Vs_from_Mu_Rho(Mu, Rho):
	"""
	Computes the shear velocity from shear modulus and density.
	The units must be as described below
	.
		Input Variables:
			Mu  = Shear Modulus          (Gpa)
			Rho = Density                (g/cc)

		Output Variable:
			Vp  = Compressional Velocity (m/s)

	"""
	# Conversion to (mks) system for computation	
	Mu = 1000000000*Mu       # Conversion from GPa to Pa
	Rho = 1000 * Rho         # Conversion from g/(cm^3) to Kg/(m^3)
	
	Vs = sqrt( Mu/Rho )
	

	return Vp
	
def get_K_from_Vp_Vs_Rho( Vp, Vs, Rho ):
	"""
	Computes bulk modulus from compressional velocity, shear velocity and density.
	The units must be as described below.
	
		Input Variables:
			Vp	= Compressional Velocity)		(m/s)
			Vs	= (Shear Velocity)				(m/s)
			Rho	= (Density)						(g/(cm^3))
			
		Output Variable:
			K	= Bulk Modulus					(Gpa)
	"""
	
	# Conversion to (mks) system for computation	
	Rho = 1000 * Rho         # Conversion from g/(cm^3) to Kg/(m^3)
	
	A = 4*Vs*Vs/3
	B = Vp*Vp
	C = B - A
	K = Rho*C
	
	K = K/1000000000         # Conversion from Pa to Gpa
	
	return K
	
def get_Mu_from_Vs_Rho( Vs, Rho ):
	"""
	Computes shear modulus from shear velocity and density.
	The units must be as described below.
		
		Input Variables:
			Vs	= (Shear Velocity)				(m/s)
			Rho	= (Density)						(g/(cm^3))
			
		Output Variable:
			Mu	= Shear Modulus					(Gpa)
	"""
	# Conversion to (mks) system for computation	
	Rho = 1000 * Rho         # Conversion from g/(cm^3) to Kg/(m^3)
	
	Mu = Rho*Vs*Vs

	Mu = Mu/1000000000       # Conversion from Pa to Gpa
	
	return K
	
def get_Ksat_FromGassmannEquation(Kdry, Kfluid, Kmin, Phi):
	"""
	Computes the saturated bulk modulus of a rock usin Gassmann's equation.
	
		Input Variables:
			Kdry = The dry frame bulk modulus			(Gpa)
			Kfluid = The fluid bulk modulus				(Gpa)
			Kmin   = The mineral/matrix bulk modulus	(Gpa)
			Phi    = The porosity 						(Fraction, V/V)
			
		Output Variable:
			Ksat   = The saturated bulk modulus			(Gpa)
			
	"""
	A1 = 1 - (Kdry/Kmin)
	A = A1 * A1
	B1 = Phi/Kfluid
	B2 = (1 - Phi)/Kmin
	B3 = Kdry/(Kmin*Kmin)
	B = B1 + B2 - B3
	C = A/B
	
	Ksat = Kfluid + C
		
	return Ksat
	
def get_Kdry_FromGassmannEquation(Ksat, Kfluid, Kmin, Phi):
	"""
	Computes the saturated bulk modulus of a rock usin Gassmann's equation.
	
		Input Variables:
			Ksat   = The saturated bulk modulus			(Gpa)
			Kfluid = The fluid bulk modulus				(Gpa)
			Kmin   = The mineral/matrix bulk modulus	(Gpa)
			Phi    = The porosity 						(Fraction, V/V)
			
		Output Variable:
			Kdry   = The dry frame bulk modulus			(Gpa)
			
	"""
	
	A1 = B1 = (Phi * Kmin)/Kfluid
	A2 = A1 + 1 - Phi
	A3 = Ksat * A2
	A = A3 - Kmin
	B2 = Ksat/Kmin
	B = B1 + B2 - 1 - Phi
	
	Kdry = A/B
	
	return Kdry
	
def voigtAverage(values=[], weights=[]):
	"""
	Computes the Voigt average (voigtAve = w1*v1 + w2*v2 + ...)
	Weights will be normalized so sum will equal 1.
	"""
	voigtAve = 0
	sumWeight = 0
	if len(values) != len(weights):
		mm.myStdErrorMessage( "Length missmatch", "The number of elements in your values does not equal the number of weights" )
		mm.exitNow()
		
	for i in len(weights):
		sumWeight += weights[i]
		
	for i in len(values):
		voigtAve += weights[i]*values[i]/sumWeight
		
	return voigtAve
	
def reussAverage(values=[], weights=[]):
	"""
	Computes the Voigt average (1/voigtAve = w1/v1 + w2/v2 + ...)
	Weights will be normalized so sum will equal 1.
	"""
	reussAve = 0
	sumWeight = 0
	if len(values) != len(weights):
		mm.myStdErrorMessage( "Length missmatch", "The number of elements in your values does not equal the number of weights" )
		mm.exitNow()
		
	for i in len(weights):
		sumWeight += weights[i]
		
	for i in len(values):
		reussAve += weights[i]/(sumWeight*values[i])
	reussAve = 1/reussAve
	
	return reussAve
	
def vrhAverage(values=[], weights=[]):
	"""
	Computes the Voigt-Reuss-Hill average (vrhAve = 1/2(voigtAve + reussAve))
	Weights will be normalized so sum will equal 1.
	"""
	reussAve = 0
	sumWeight = 0
	if len(values) != len(weights):
		mm.myStdErrorMessage( "Length missmatch", "The number of elements in your values does not equal the number of weights" )
		mm.exitNow()
		
	voigtAve = voigtAverage( values, weights )
	reussAve = reussAverage( values, weights )
	vrhave   = 0.5*(voigtAve + reussAve)
		
	return vrhAve
	
def hashinShtrikmanK( K1, Mu1, f1, K2, Mu2, f2 ):
	"""
	Computes Hashin Shtrikman bounds on the bulk modulus for 2 material mixture.
	Upper bounds is acheived when the subscrpt 1 is populated with the stiff 
	material and visa versa. Definition of variables follows.
		
		Input Variables:
			K1	=	Bulk modulus of one phase of the mixture		(GPa)
			Mu1	=	Shear modulus of one phase of the mixture		(GPa)
			f1	=	Volume fraction of one phase of the mixture		(fraction, V/V)
			K2	=	Bulk modulus of the other phase of the mixture	(GPa)
			Mu2	=	Shear modulus of the other phase of the mixture	(GPa)
			f1	=	Volume fraction of one phase of the mixture		(fraction, V/V)
		
		Output Variable:
			Khs	=	Hashin Shtrikman bulk modulus					(GPa)
			
	Recall if (K1,Mu1,f1) corresponds to the stiffer material this will be an upper bound,
	otherwise, if (K1,Mu1,f1) corresponds to the softer material this will be the lower bound.
	
	Mavko-Stanford Rock Physics Lab
	"""
	A = 1/(K2-K1)
	B = f1/(K1+3*Mu1/3)
	C = f2/(A+B)
	Khs = K1 + C
	
	return Khs
	
def hashinShtrikmanMu( K1, Mu1, f1, K2, Mu2, f2 ):
	"""
	Computes Hashin Shtrikman bounds on the shear modulus for 2 material mixture.
	Upper bounds is acheived when the subscrpt 1 is populated with the stiff 
	material and visa versa. Definition of variables follows.
		
		Input Variables:
			K1		=	Bulk modulus of first phase of the mixture		(GPa)
			Mu1		=	Shear modulus of first phase of the mixture		(GPa)
			f1		=	Volume fraction of first phase of the mixture	(fraction, V/V)
			K2		=	Bulk modulus of the 2nd phase of the mixture	(GPa)
			Mu2		=	Shear modulus of the 2nd phase of the mixture	(GPa)
			f1		=	Volume fraction of second phase of the mixture	(fraction, V/V)
		
		Output Variable:
			Muhs	=	Hashin Shtrikman bulk modulus					(GPa)
			
	Recall if (K1,Mu1,f1) corresponds to the stiffer material this will be an upper bound,
	otherwise, if (K1,Mu1,f1) corresponds to the softer material this will be the lower bound.
	
	Mavko-Stanford Rock Physics Lab
	"""
	A = 1/(Mu2-Mu1)
	B1 = 2*f1*(K1+2*Mu1)
	B2 = 5*Mu1*(K1+4*Mu1/3)
	B = B1/B2
	C = f2/(A+B)
	
	Muhs = Mu1 + C
	
	return Muhs
	
def applyScalarFunctionToArrayVariables( function, parms=[]):
	"""
	apply a scalar function to an array of variables each of which can be arrays as well
	right now this function can only handle up to 9 variable arrays of equal dimention.
	
	"""
	numvars = len(parms)
	
	if numvars == 1:
		output = [function(parms[0][i]) for i in range(len(parms[0]))]
	elif numvars == 2:
		output = [function(parms[0][i],parms[1][i]) for i in range(len(parms[0]))]
	elif numvars == 3:
		output = [function(parms[0][i],parms[1][i],parms[2][i]) for i in range(len(parms[0]))]
	elif numvars == 4:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i]) for i in range(len(parms[0]))]
	elif numvars == 5:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i],parms[4][i]) for i in range(len(parms[0]))]
	elif numvars == 6:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i],parms[4][i],parms[5][i]) for i in range(len(parms[0]))]
	elif numvars == 7:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i],parms[4][i],parms[5][i],parms[6][i]) for i in range(len(parms[0]))]
	elif numvars == 8:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i],parms[4][i],parms[5][i],parms[6][i],parms[7][i]) for i in range(len(parms[0]))]
	elif numvars == 9:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i],parms[4][i],parms[5][i],parms[6][i],parms[7][i],parms[8][i]) for i in range(len(parms[0]))]
	elif numvars == 10:
		output = [function(parms[0][i],parms[1][i],parms[2][i],parms[3][i],parms[4][i],parms[5][i],parms[6][i],parms[7][i],parms[8][i]) for i in range(len(parms[0]))]
	else:
		mm.myStdErrorMessage( "Too many variables", "Only up to 10 variables can be entered for a funtion" )
		mm.exitNow()
		
	return output

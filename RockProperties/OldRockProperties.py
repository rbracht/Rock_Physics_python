# My Python rock properties modules

import numpy as np
import scipy as sp
import MyModule as mm
from math import sqrt, exp

def get_Vp_from_K_Mu_Rho(K, Mu, Rho ):
	"""
	Computes the compressional velocity from bulk modulus, shear modulus and density.
	The units must be as described below.
	
		Input Variables:
			K   = Bulk modulus           (Gpa)
			Mu  = Shear Modulus          (Gpa)
			Rho = Density                (g/cc)
			
		Output Variable:
			Vp  = Compressional Velocity (m/s)
			
	Note the inputs variables can be either scalars, python lists or numpy arrays, 
	however all variables must be of the same type and dimenson in the case of lists
	and arrays.
	"""
	if isinstance(K,(int, long, float, complex)) and (
		 isinstance(Mu,(int, long, float, complex))) and (
		 isinstance(Rho,(int, long, float, complex))):
					
		# Conversion to (mks) system for computation
		K = 1E9 * K         # Conversion from GPa to Pa
		Mu = 1E9 * Mu       # Conversion from GPa to Pa
		Rho = 1E3 * Rho     # Conversion from g/(cm^3) to Kg/(m^3)
		
		Vp = sqrt( (K + 4*Mu/3)/Rho )
		
	elif type(K) == np.ndarray and (
		type(Mu) == np.ndarray) and (
		type(Rho) == np.ndarray) and (
		len(K) == len(Mu)) and (
		len(Mu) == len(Rho)):
		
		Vp = np.sqrt( (1E9*K + 4E9*Mu/3)/(1E3*Rho) )
		
	elif type(K) == list and (
		type(Mu) == list) and (
		type(Rho) == list) and (
		len(K) == len(Mu)) and (
		len(Mu) == len(Rho)):
		
		Vp = [sqrt((1E9*K[i] + 4E9*Mu[i]/3)/(1E3*Rho[i])) for i in range(len(K))]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", 
		"The input variables are not consistent for this function" )
		mm.exitNow()
	
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
			
	Note the inputs variables can be either scalars, python lists or numpy arrays, 
	however all variables must be of the same type and dimenson in the case of lists
	and arrays.
	"""
	if isinstance(Mu,(int, long, float, complex)) and (
		isinstance(Rho,(int, long, float, complex))):
		
		# Conversion to (mks) system for computation
		Mu = 1E9 * Mu       # Conversion from GPa to Pa
		Rho = 1E3 * Rho     # Conversion from g/(cm^3) to Kg/(m^3)
		
		Vs = sqrt( Mu/Rho )
		
	elif type(Mu) == np.ndarray and (
		type(Rho) == np.ndarray) and (
		len(Mu) == len(Rho)):

		Vs = np.sqrt( (1E9*Mu)/(1E3*Rho) )
		
	elif type(Mu) == list and (
		type(Rho) == list) and (
		len(Mu) == len(Rho)):
		
		Vs = [sqrt((1E9*Mu[i])/(1E3*Rho[i])) for i in range(len(Mu))]
		
	else:

		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
		
	return Vs
	
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
	
	if isinstance(Vp,(int, long, float, complex)) and (
		isinstance(Vs,(int, long, float, complex))) and (
		isinstance(Rho,(int, long, float, complex))):
		
		# Conversion to (mks) system for computation	
		Rho = 1E3 * Rho   # Conversion from g/(cm^3) to Kg/(m^3)

		A = 4*Vs*Vs/3
		B = Vp*Vp
		C = B - A
		K = Rho*C

		K = K/1E9         # Conversion from Pa to Gpa
		
	elif type(Vp) == np.ndarray and (
		type(Vs) == np.ndarray) and (
		type(Rho) == np.ndarray) and (
		len(Vp) == len(Vs)) and (
		len(Vs) == len(Rho)):
		
		# Conversion to (mks) system for computation	
		Rho = 1E3 * Rho    # Conversion from g/(cm^3) to Kg/(m^3)

		A = 4*(Vs*Vs)/3
		B = Vp*Vp
		C = B - A
		K = Rho*C

		K = K/1E9         # Conversion from Pa to Gpa
		
	elif type(Vp) == list and (
		type(Vs) == list) and (
		type(Rho) == list) and (
		len(Vp) == len(Vs)) and (
		len(Vs) == len(Rho)):
		
		Rho = [Rho[i]*1E3 for i in range(len(Rho))] # Conversion from g/(cm^3) to Kg/(m^3)
		
		A = [4*Vs[i]*Vs[i]/3 for i in range(len(Vs))]
		B = [Vp[i]*Vp[i] for i in range(len(Vp))]
		C = [B[i]-A[i] for i in range(len(B))]
		K = [Rho[i] * C[i] for i in range(len(C))]
		
		K = [K[i]/1E9 for i in range(len(K))] # Conversion from Pa to Gpa
		
	else:

		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
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
	
	if isinstance(Vs,(int, long, float, complex)) and (
		isinstance(Rho,(int, long, float, complex))):
		
		# Conversion to (mks) system for computation	
		Rho = 1E3 * Rho   # Conversion from g/(cm^3) to Kg/(m^3)
		
		Mu = Rho*(Vs*Vs)
		
		Mu = Mu/1E9         # Conversion from Pa to Gpa
		
	elif type(Vs) == np.ndarray and (
		type(Rho) == np.ndarray) and (
		len(Vs) == len(Rho)):
		
		# Conversion to (mks) system for computation	
		Rho = 1E3 * Rho    # Conversion from g/(cm^3) to Kg/(m^3)
		
		Mu = Rho*(Vs*Vs)
		
		Mu = Mu/1E9         # Conversion from Pa to Gpa
		
	elif type(Vs) == list and (
		type(Rho) == list) and (
		len(Vs) == len(Rho)):
		
		Rho = [Rho[i]*1E3 for i in range(len(Rho))] # Conversion from g/(cm^3) to Kg/(m^3)
		
		Mu = [Rho[i] * (Vs[i]*Vs[i]) for i in range(len(Vs))]
		
		Mu = [Mu[i]/1E9 for i in range(len(Mu))] # Conversion from Pa to Gpa
		
	else:

		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
	return Mu
	
def get_Ksat_FromGassmannEquation(Kdry, Kfluid, Kmin, Phi):
	"""
	Computes the saturated bulk modulus of a rock usin Gassmann's equation.
	
		Input Variables:
			Kdry	= The dry frame bulk modulus			(Gpa)
			Kfluid	= The fluid bulk modulus				(Gpa)
			Kmin	= The mineral/matrix bulk modulus		(Gpa)
			Phi		= The porosity 							(Fraction, V/V)
			
		Output Variable:
			Ksat	= The saturated bulk modulus			(Gpa)
			
	"""
	
	if isinstance(Kdry,(int, long, float, complex)) and (
		isinstance(Kfluid,(int, long, float, complex))) and (
		isinstance(Kmin,(int, long, float, complex))) and (
		isinstance(Phi,(int, long, float, complex))):
		
		A1 = 1 - (Kdry/Kmin)
		A = A1 * A1
		B1 = Phi/Kfluid
		B2 = (1 - Phi)/Kmin
		B3 = Kdry/(Kmin*Kmin)
		B = B1 + B2 - B3
		C = A/B

		Ksat = Kfluid + C
		
	elif type(Kdry) == np.ndarray and (
		type(Kfluid) == np.ndarray) and (
		type(Kmin) == np.ndarray) and (
		type(Phi) == np.ndarray) and (
		len(Kdry) == len(Kfluid)) and (
		len(Kfluid) == len(Kmin)) and (
		len(Kmin) == len(Phi)):
		
		A1 = 1 - (Kdry/Kmin)
		A = A1 * A1
		B1 = Phi/Kfluid
		B2 = (1 - Phi)/Kmin
		B3 = Kdry/(Kmin*Kmin)
		B = B1 + B2 - B3
		C = A/B
		
		Ksat = Kfluid + C
		
	elif type(Vp) == list and (
		type(Vs) == list) and (
		type(Rho) == list) and (
		len(Vp) == len(Vs)) and (
		len(Vs) == len(Rho)):
		
		numItems = len(Vp)
		
		A1 = [1-(Kdry[i]/Kmin[i]) for i in range(numItems)]
		A = [A1[i]*A1[i] for i in range(numItems)]
		B1 = [Phi[i]/Kfluid[i] for i in range(numItems)]
		B2 = [(1-Phi[i])/Kmin[i] for i in range(numItems)]
		B3 = [Kdry[i]/(Kmin[i]*Kmin[i]) for i in range(numItems)]
		B = [B1[i] + B[i] - B3[i] for i in range(numItems)]
		C = [A[i]/B[i] for i in range(numItems)]
		
		Ksat = [Kfluid[i] + C[i] for i in range(numItems)]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
	return Ksat
	
def get_Kdry_FromGassmannEquation(Ksat, Kfluid, Kmin, Phi):
	"""
	Computes the saturated dry frame bulk modulus of a rock using Gassmann's equation.
	
		Input Variables:
			Ksat   = The saturated bulk modulus			(Gpa)
			Kfluid = The fluid bulk modulus				(Gpa)
			Kmin   = The mineral/matrix bulk modulus	(Gpa)
			Phi    = The porosity 						(Fraction, V/V)
			
		Output Variable:
			Kdry   = The dry frame bulk modulus			(Gpa)
			
	"""
	
	if isinstance(Ksat,(int, long, float, complex)) and (
		isinstance(Kfluid,(int, long, float, complex))) and (
		isinstance(Kmin,(int, long, float, complex))) and (
		isinstance(Phi,(int, long, float, complex))):
		
		A1 = B1 = (Phi * Kmin)/Kfluid
		A2 = A1 + 1 - Phi
		A3 = Ksat * A2
		A = A3 - Kmin
		B2 = Ksat/Kmin
		B = B1 + B2 - 1 - Phi

		Kdry = A/B
		
	elif type(Kdry) == np.ndarray and (
		type(Kfluid) == np.ndarray) and (
		type(Kmin) == np.ndarray) and (
		type(Phi) == np.ndarray) and (
		len(Kdry) == len(Kfluid)) and (
		len(Kfluid) == len(Kmin)) and (
		len(Kmin) == len(Phi)):
		
		A1 = B1 = (Phi * Kmin)/Kfluid
		A2 = A1 + 1 - Phi
		A3 = Ksat * A2
		A = A3 - Kmin
		B2 = Ksat/Kmin
		B = B1 + B2 - 1 - Phi

		Kdry = A/B
		
	elif type(Vp) == list and (
		type(Vs) == list) and (
		type(Rho) == list) and (
		len(Vp) == len(Vs)) and (
		len(Vs) == len(Rho)):
		
		numItems = len(Vp)
		
		A1 = B1 = [(Phi[i] * Kmin[i]) / Kfluid[i] for i in range(numItems)]
		A2 = [A1[i] + 1 - Phi[i] for i in range(numItems)]
		A3 = [Ksat[i] * A2[i] for i in range(numItems)]
		A = [A3[i] - Kmin[i] for i in range(numItems)]
		B2 = [Ksat[i] / Kmin[i] for i in range(numItems)]
		B = [B1[i] + B2[i] - 1 - Phi[i] for i in range(numItems)]
		
		Kdry = [A[i] / B[i] for i in range(numItems)]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
	Kdry = A/B
	
	return Kdry
	
def voigtAverage(values=[], weights=[]):
	"""
	Computes the Voigt average (voigtAve = w1*v1 + w2*v2 + ...)
	Weights will be normalized so sum will equal 1.
	"""
	voigtAve = 0
	sumWeight = 0
	i = 0
	
	if type(values) == np.ndarray and (
		type(weights) == np.ndarray) and (
		len(values) == len(weights)):
		
		for w in np.nditer(weights):
			sumWeight += w
			
		for i in range(len(values)):
			voigtAve += weights[i]*values[i]/sumWeight
		
	elif type(values) == list and (
		type(weights) == list) and (
		len(values) == len(weights)):
		
		for i in range(len(weights)):
			sumWeight += weights[i]
			
		for i in range(len(values)):
			voigtAve += weights[i]*values[i]/sumWeight
		
	else:
		
		mm.myStdErrorMessage( "Error: voigtAverage", (
			"The input variables are not consistent for this function" ))
		return
	
	return voigtAve
	
def reussAverage(values=[], weights=[]):
	"""
	Computes the Voigt average (1/reussAve = w1/v1 + w2/v2 + ...)
	Weights will be normalized to insure the sum of all weights will equal 1.
	It is assumed, both input values are 1 dimensional.
	"""
	reussAve = 0
	sumWeight = 0
	i = 0
	
	if type(values) == np.ndarray and (
		type(weights) == np.ndarray) and (
		len(values) == len(weights)):
		
		for i in range(len(weights)):
			sumWeight += weights[i]

		for i in range(len(values)):
			reussAve += weights[i]/(values[i]*sumWeight)
		
		reussAve = 1/reussAve
		
	elif type(values) == list and (
		type(weights) == list) and (
		len(values) == len(weights)):
		
		for i in range(len(weights)):
			sumWeight += weights[i]

		for i in range(len(values)):
			reussAve += weights[i]/(values[i]*sumWeight)
			
		reussAve = 1/reussAve
		
	else:
		
		mm.myStdErrorMessage( "Error: voigtAverage", (
			"The input variables are not consistent for this function" ))
		return
	
	return reussAve
	
def vrhAverage(values=[], weights=[]):
	"""
	Computes the Voigt-Reuss-Hill average (vrhAve = 1/2(voigtAve + reussAve))
	Weights will be normalized so sum will equal 1.
	"""
	reussAve = 0
	sumWeight = 0
	if len(values) != len(weights):
		mm.myStdErrorMessage( "Error: vrhAverage", (
		 "The number of elements in your values does not equal the number of weights" ))
		mm.exitNow()
		
	voigtAve = voigtAverage( values, weights )
	reussAve = reussAverage( values, weights )
	vrhAve   = 0.5*(voigtAve + reussAve)
		
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
	
	if isinstance(K1,(int, long, float, complex)) and (
		isinstance(Mu1,(int, long, float, complex))) and (
		isinstance(f1,(int, long, float, complex))) and (
		isinstance(K2,(int, long, float, complex))) and (
		isinstance(Mu2,(int, long, float, complex))) and (
		isinstance(f2,(int, long, float, complex))):
		
		A = 1/(K2-K1)
		B = f1/(K1+4*Mu1/3)
		C = f2/(A+B)
		Khs = K1 + C
		
	elif type(K1) == np.ndarray and (
		type(Mu1) == np.ndarray) and (
		type(f1) == np.ndarray) and (
		type(K2) == np.ndarray) and (
		type(Mu2) == np.ndarray) and (
		type(f2) == np.ndarray) and (
		len(K1) == len(Mu1)) and (
		len(Mu1) == len(f1)) and (
		len(f1) == len(K2)) and (
		len(K2) == len(Mu2)) and (
		len(Mu2) == len(f2)):
		
		A = 1/(K2-K1)
		B = f1/(K1+4*Mu1/3)
		C = f2/(A+B)
		Khs = K1 + C
		
	elif type(K1) == list and (
		type(Mu1) == list) and (
		type(f1) == list) and (
		type(K2) == list) and (
		type(f2) == list) and (
		len(K1) == len(Mu1)) and (
		len(Mu1) == len(f1)) and (
		len(f1) == len(K2)) and (
		len(K2) == len(Mu2)) and (
		len(Mu2) == len(f2)):
		
		numItems = len(K1)
		
		A = [1/(K2[i] * K1[i]) for i in range(numItems)]
		B = [f1[i]/(K1[i] + (4*Mu1[i]/3)) for i in range(numItems)]
		C = [f2[i] / (A[i]+B[i]) for i in range(numItems)]
		
		Khs = [K1[i] + C[i] for i in range(numItems)]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
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
	
	if isinstance(K1,(int, long, float, complex)) and (
		isinstance(Mu1,(int, long, float, complex))) and (
		isinstance(f1,(int, long, float, complex))) and (
		isinstance(K2,(int, long, float, complex))) and (
		isinstance(Mu2,(int, long, float, complex))) and (
		isinstance(f2,(int, long, float, complex))):
		
		A = 1/(Mu2-Mu1)
		B1 = 2*f1*(K1+2*Mu1)
		B2 = 5*Mu1*(K1+4*Mu1/3)
		B = B1/B2
		C = f2/(A+B)

		Muhs = Mu1 + C
		
	elif type(K1) == np.ndarray and (
		type(Mu1) == np.ndarray) and (
		type(f1) == np.ndarray) and (
		type(K2) == np.ndarray) and (
		type(Mu2) == np.ndarray) and (
		type(f2) == np.ndarray) and (
		len(K1) == len(Mu1)) and (
		len(Mu1) == len(f1)) and (
		len(f1) == len(K2)) and (
		len(K2) == len(Mu2)) and (
		len(Mu2) == len(f2)):
		
		A = 1/(Mu2-Mu1)
		B1 = 2*f1*(K1+2*Mu1)
		B2 = 5*Mu1*(K1+4*Mu1/3)
		B = B1/B2
		C = f2/(A+B)
		
		Muhs = Mu1 + C
		
	elif type(K1) == list and (
		type(Mu1) == list) and (
		type(f1) == list) and (
		type(K2) == list) and (
		type(f2) == list) and (
		len(K1) == len(Mu1)) and (
		len(Mu1) == len(f1)) and (
		len(f1) == len(K2)) and (
		len(K2) == len(Mu2)) and (
		len(Mu2) == len(f2)):
		
		numItems = len(K1)
		
		A = [1/(Mu2[i] - Mu1[i]) for i in range(numItems)]
		B1 = [2*f1[i]*(K1[i] + (2*Mu1[i])) for i in range(numItems)]
		B2 = [5*Mu1[i]*(K1[i] + (4*Mu1[i]/3)) for i in range(numItems)]
		B = [B1[i] / B2[i] for i in range(numItems)]
		C = [f2[i] / (A[i]+B[i]) for i in range(numItems)]
		
		Muhs = [Mu1[i] + C[i] for i in range(numItems)]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
	return Muhs
	
def applyScalarFunctionToArrayVariables( function, parms=[]):
	"""
	apply a scalar function to an array of variables each of which can be arrays as well
	right now this function can only handle up to 10 variable arrays of equal dimension.
	
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
	
def matrix_K_calculator(constituentVolFrac = [], constituentK = []):
	"""
	The function computes the equivalent bulk modulus Kmin from the 
	contituent minerals within the matrix and their respective volume fractions
	using the Voigt-Reuss-Hill average.
	"""
	
	if len(constituentVolFrac) == len(constituentK):
	
		Kmin = vrhAverage(constituentVolFrac, constituentK)
		
	else:
		mm.myStdErrorMessage( "Length missmatch", 
		"The number of elements in weights do not match the numvber of constituents" )
		return
		
	return Kmin
	
def brineDensity( salinity, reservoirP, reservoirT ):
	"""
	This function computes the brine density at reservoir conditions.
	The only input required are:
	
		salinity:		the salinity of formation brine (weight fraction)
		reservoirP:		reservoir pressure (MPa)
		reservoirT:		reservoir temperature (degrees C)
		
	The output will be:
	
		brineDensity:	formation brine density (g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	S = salinity
	
	A1 = -80*T-3.3*T*T+0.00175*T*T*T+489*P-2*T*P
	A2 = 0.016*T*T*P-1.3E-5*T*T*T*P-0.333*P*P-0.002*T*P*P
	waterDensity = 1 +1E-6*(A1 + A2)
	
	B1 = 80+3*T-3300*S-13*P+47*P*S
	B2 = 300*P-2400*P*S
	B3 = 0.668*S + 0.44*S*S

	brineDensity = waterDensity + B3 + 1E-6*S*(B2 + T*B1)

	return brineDensity
	
def pVelBrine(salinity, reservoirP, reservoirT):
	"""
	The function computes the compressional velocity of brine at reservoir conditions.
	The only input required are:
	
		salinity:		the salinity of formation brine (weight fraction)
		reservoirP:		reservoir pressure (MPa)
		reservoirT:		reservoir temperature (degrees C)
		
	The output will be:
	
		pVelBrine:		formation brine velocity (m/s)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	S = salinity

	w = [[0 for i in range(5)] for j in range(6)]
	
	w[1][1] = 1402.85
	w[2][1] = 4.871
	w[3][1] = -0.04783
	w[4][1] = 1.487E-4
	w[5][1] = -2.197E-7
	w[1][2] = 1.524
	w[2][2] = -0.0111
	w[3][2] = 2.747E-4
	w[4][2] = -6.503E-7
	w[5][2] = 7.987E-10
	w[1][3] = 3.437E-3
	w[2][3] = 1.739E-4
	w[3][3] = -2.1357E-6
	w[4][3] = -1.455E-8
	w[5][3] = 5.23E-11
	w[1][4] = -1.197E-5
	w[2][4] = -1.628E-6
	w[3][4] = 1.237E-8
	w[4][4] = 1.327E-10
	w[5][4] = -4.614E-13
	
	pVelWater = 0.0
	for i in range(1,6):
		for j in range(1,5):
			pVelWater += w[i][j]*T**(i-1)*P**(j-1)
			
	A = 1170-9.6*T+0.055*T*T-8.5E-5*T*T*T+2.6*P-0.0029*T*P-0.0476*P*P
	B = 780-10*P+0.16*P*P
	pVelBrine = pVelWater + S*A + S**1.5*B - 1820*S*S
	
	return pVelBrine
	
	
def brineBulkModulus(salinity, reservoirP, reservoirT):
	"""
	This function computes the bulk mudulus of brine at reservoir conditions.
	The only input required are:
	
		salinity:		the salinity of formation brine (weight fraction)
		reservoirP:		reservoir pressure (MPa)
		reservoirT:		reservoir temperature (degrees C)
		
	The output will be:
	
		Kbrine:			Brine bulk modulus at reservoir conditions (GPa)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	S = salinity
	
	Rbrine = brineDensity(S,P,T)
	Vbrine = pVelBrine(S,P,T)
	
	Kbrine = Rbrine*Vbrine*Vbrine*1E-6
	
	return Kbrine
	
def gasDensity( gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the density of gas at reservoir conditions.
	The only input required are:
	
		gasSpecificGravity:      	the specific gravity of gas (API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure (MPa)
		reservoirT:			reservoir temperature (degrees C)
		
	The output will be:
	
		Rgas:				gas density at reservoir conditions (g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	R = 8.314  #Gas constant
	
	Ppr = P/(4.892 - 0.4048*G)
	Tpr = (T+273.15)/(94.72 + 170.75*G)
	
	A = (0.45+8*(0.56-1/Tpr)**2)*(-Ppr**1.2/Tpr)
	B = 0.109*(3.85-Tpr)**2
	E = B*exp(A)
	
	C = 0.03 + 0.00527*(3.5 - Tpr)**3
	D = 0.642*Tpr-0.007*Tpr**4-0.52
	Z = C*Ppr + D + E
	
	Rgas = (28.8*G*P)/(Z*R*(T+273.15))
	
	return Rgas

def gasBulkModulus( gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the density of gas at reservoir conditions.
	The only input required are:
	
		gasSpecificGravity:      	the specific gravity of gas (API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure (MPa)
		reservoirT:			reservoir temperature (degrees C)
		
	The output will be:
	
		Rgas:				gas density at reservoir conditions (g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	R = 8.314  #Gas constant
	
	Ppr = P/(4.892 - 0.4048*G)
	Tpr = (T+273.15)/(94.72 + 170.75*G)
	
	A = (0.45+8*(0.56-1/Tpr)**2)*(-Ppr**1.2/Tpr)
	B = 0.109*(3.85-Tpr)**2
	E = B*exp(A)
	
	C = 0.03 + 0.00527*(3.5 - Tpr)**3
	D = 0.642*Tpr-0.007*Tpr**4-0.52
	Z = C*Ppr + D + E

	A1 = -1.2*(Ppr**0.2/Tpr)*(0.45+8*(0.56-1/Tpr)**2)
	B1 = (0.45+8*(0.56-1/Tpr)**2)*(-Ppr**1.2/Tpr)
	F = A1*exp(B1)
	
	dZdP = 0.03 + 0.00527*(3.5-Tpr)**3+0.109*(3.85-Tpr)**2 * F
	Gama0 = 0.85+(5.6/(Ppr+2))+(27.1/(Ppr+3.5)**2)-8.7*exp(-0.65*(Ppr+1))

	Kgas = (Gama0/1000)*P/(1-Ppr*dZdP/Z)

	return Kgas

def oilDensity( oilApiGravity, gasOilRatio, gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the oil denisty at reservoir conditions.
	The only input required are:

		oilApiGravity:          	API gravity of oil  (API at 15.6C and 1 atmosphere)
	        gasOilRatio:                    gas to oil ratio (litre/litre)
		gasSpecificGravity:      	ratio of gas density to air density (API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure (MPa)
		reservoirT:			reservoir temperature (degrees C)
		
	The output will be:
	
		Roil:				oil density at reservoir conditions (g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	Rg = gasOilRatio
	R0 = 141.5/( oilApiGravity + 131.5 )     # Convert from api gravity to specific gravity gm/cc

	B1 = (2.495 * Rg * sqrt(G/R0) + T + 17.8)**1.175 
	B0 = 0.972 + 0.00038 * B1                # formation volume factor

	Rps = R0 / ((1 + 0.001*Rg) * B0)
	Rs = (R0 + 0.0012*Rg*G) / B0

	C1 = (0.00277*P - 1.71E-7 * P*P*P) * (Rs - 1.15)**2
	C2 = 3.49E-4 * P
	C3 = Rs + C1 + C2
	C4 = 0.972 + 3.81E-4*(T + 17.78)**1.175
	
	Roil = C3 / C4

	return Roil

def oilPwaveVelocity( oilApiGravity, gasOilRatio, gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the oil denisty at reservoir conditions.
	The only input required are:

		oilApiGravity:          	API gravity of oil  (API at 15.6C and 1 atmosphere)
	        gasOilRatio:                    gas to oil ratio (litre/litre)
		gasSpecificGravity:      	ratio of gas density to air density (API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure (MPa)
		reservoirT:			reservoir temperature (degrees C)
		
	The output will be:
	
		Vp_oil:				oil p-wave velocity at reservoir conditions (g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	Rg = gasOilRatio
	R0 = 141.5/( oilApiGravity + 131.5 )   # Convert from api gravity to specific gravity gm/cc

	B1 = (2.495 * Rg * sqrt(G/R0) + T + 17.8)**1.175
	B0 = 0.972 + 0.00038 * B1                  # formation volume factor

	Rps = R0 / ((1 + 0.001 * Rg) * B0)
	Rs = (R0 + 0.0012 * Rg * G) / B0

	A = 2096*sqrt(Rps / (2.6 - Rps)) - 3.7*T + 4.64*P
	B = 0.0115*(sqrt((18.33/Rps) - 16.97) - 1)*T*P
	Vp_oil = A + B

	return Vp_oil

def oilBulkModulus(oilApiGravity, gasOilRatio, gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the oil denisty at reservoir conditions.
	The only input required are:

		oilApiGravity:          	API gravity of oil  (API at 15.6C and 1 atmosphere)
	        gasOilRatio:                    gas to oil ratio (litre/litre)
		gasSpecificGravity:      	ratio of gas density to air density (API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure (MPa)
		reservoirT:			reservoir temperature (degrees C)
		
	The output will be:
	
		Koil:				oil bulk modulus at reservoir conditions (Gpa)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	Rg = gasOilRatio
	R0 = 141.5/( oilApiGravity + 131.5 )    # convert from api to gm/cc

	Voil = oilPwaveVelocity(R0,Rg,G,P,T)    # Oil p-wave velocity m/s
	Roil = oilDensity(R0,Rg,G,P,T)          # Oll density g/cm^3

	Koil = Roil*Voil*Voil*1E-6

	return Koil

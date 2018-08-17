# My Python rock properties modules
# Updated 8 September 2013 Don Easley on Mac at home

import numpy as np
import scipy as sp
from scipy.optimize import fsolve, bisect
import MyModule as mm
import matplotlib.pyplot as plt
from math import sqrt, exp, pi, fabs, log

def getPoissonRatioFromBulkAndShearModulus( bulkModulus, shearModulus ):
	"""
	The function computes the Poisson's ratio (Nu) given a bulk modulus (K)
	and a shear modulus (Mu). Make sure the units for both bulk and shear 
	modulus are the same (such as GPa). The Poisson's ratio is unitless.
	"""
	K = bulkModulus
	Mu = shearModulus
	
	Nu = (3.0*K - 2.0*Mu)/(2.0*(3*K + Mu))
	
	return Nu

def get_Vp_from_K_Mu_Rho( K, Mu, Rho ):
	"""
	Computes the compressional velocity from bulk modulus, shear modulus and density.
	The units must be as described below.
	
		Input Variables:
			K   = Bulk modulus      	(Gpa)
			Mu  = Shear Modulus      	(Gpa)
			Rho = Density		        (gm/cc)
			
		Output Variable:
			Vp  = Compressional Velocity	(m/s)
			
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
	
def get_Vs_from_Mu_Rho( Mu, Rho ):
	"""
	Computes the shear velocity from shear modulus and density.
	The units must be as described below
	.
		Input Variables:
			Mu  = Shear Modulus		(Gpa)
			Rho = Density			(g/cc)

		Output Variable:
			Vp  = Compressional Velocity	(m/s)
			
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
			Vp	= Compressional Velocity)	(m/s)
			Vs	= (Shear Velocity)		(m/s)
			Rho	= (Density)			(g/(cm^3))
			
		Output Variable:
			K	= Bulk Modulus			(Gpa)
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
			Vs	= (Shear Velocity)	(m/s)
			Rho	= (Density)		(g/(cm^3))
			
		Output Variable:
			Mu	= Shear Modulus		(Gpa)
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
	
def get_Ksat_FromGassmannEquation( Kdry, Kfluid, Kmin, Phi ):
	"""
	Computes the saturated bulk modulus of a rock usin Gassmann's equation.
	
		Input Variables:
			Kdry	= The dry frame bulk modulus		(Gpa)
			Kfluid	= The fluid bulk modulus		(Gpa)
			Kmin	= The mineral/matrix bulk modulus	(Gpa)
			Phi		= The porosity 			(Fraction, V/V)
			
		Output Variable:
			Ksat	= The saturated bulk modulus		(Gpa)
			
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

		Ksat = Kdry + C
		
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
		
		Ksat = Kdry + C
		
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
		
		Ksat = [Kdry[i] + C[i] for i in range(numItems)]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
	return Ksat
	
def get_Kdry_FromGassmannEquation( Ksat, Kfluid, Kmin, Phi, phiCutoff=0.05, null=-999.25 ):
	"""
	Computes the saturated dry frame bulk modulus of a rock using Gassmann's equation.
	A porosity cutoff is incorporated so any porosity less than this value will result
	in the dry frame modulus being equal to Ksat
	
		Input Variables:
			Ksat   = The saturated bulk modulus		(Gpa)
			Kfluid = The fluid bulk modulus			(Gpa)
			Kmin   = The mineral/matrix bulk modulus	(Gpa)
			Phi    = The porosity 				(Fraction, V/V)
			
		Output Variable:
			Kdry   = The dry frame bulk modulus		(Gpa)
			
	"""
	
	if isinstance(Ksat,(int, long, float, complex)) and (
		isinstance(Kfluid,(int, long, float, complex))) and (
		isinstance(Kmin,(int, long, float, complex))) and (
		isinstance(Phi,(int, long, float, complex))):
		
		A1 = Ksat / (Kmin - Ksat)

		if Phi == 0.0 or Phi < phiCutoff:
			Kdry = Ksat
		else:
			A1 = Ksat / (Kmin - Ksat)
			A2 = Kfluid / (Phi * (Kmin - Kfluid))
			A  = A1 - A2
			a = A / (1 + A)
			if a < 0.0 or a > 1.0:
				Kdry = null
			else:
				Kdry = a * Kmin	
			
		return Kdry
			
		
	elif type(Ksat) == np.ndarray and (
		type(Kfluid) == np.ndarray) and (
		type(Kmin) == np.ndarray) and (
		type(Phi) == np.ndarray) and (
		len(Ksat) == len(Kfluid)) and (
		len(Kfluid) == len(Kmin)) and (
		len(Kmin) == len(Phi)):
		
		A1 = B1 = (Phi * Kmin) / Kfluid
		A2 = A1 + 1.0 - Phi
		A3 = Ksat * A2
		A = A3 - Kmin
		B2 = Ksat/Kmin
		B = B1 + B2 - 1.0 - Phi

		Kdry = A/B
			
	elif type(Ksat) == list and (
		type(Kfluid) == list) and (
		type(Kmin) == list) and (
		type(Phi) == list) and (
		len(Ksat) == len(Kfluid)) and (
		len(Kfluid) == len(Kmin)) and (
		len(Kmin) == len(Phi)):
		
		numItems = len(Phi)
		
		A1 = B1 = [(Phi[i] * Kmin[i]) / Kfluid[i] for i in range(numItems)]
		A2 = [A1[i] + 1.0 - Phi[i] for i in range(numItems)]
		A3 = [Ksat[i] * A2[i] for i in range(numItems)]
		A = [A3[i] - Kmin[i] for i in range(numItems)]
		B2 = [Ksat[i] / Kmin[i] for i in range(numItems)]
		B = [B1[i] + B2[i] - 1.0 - Phi[i] for i in range(numItems)]
		
		Kdry = [A[i] / B[i] for i in range(numItems)]
		
	else:
		
		mm.myStdErrorMessage( "Error: get_Vp_from_K_Mu_Rho", (
			"The input variables are not consistent for this function" ))
		return
	
	return Kdry

def hertzMindlinSpherePackEffectiveBulkModulus(	effectivePressure,
						aveNumContactPerSphere = 9.0, 
						porosity = 0.36, 
						sphereBulkModulus = 36.6, 
						sphereShearModulus = 45.0, 
						alpha = 1.0/3.0 ):
	"""
	The function calculated the effective bulk modulus (Keff) of a random packing of spheres
	with average number of contacts per sphere (C), inter-sphere porosity (Phi), the shear
	and bulk of the spheres (Mu and K) along with the effective pressure (P),
	
	(Variable)		(Brief description)				(Units/Values)
	aveNumContactPerSphere	average number of contacts per sphere		(between 5 and 9 is common)
	porosity		inter-sphere space as fraction of whole		(fraction V/V)
	sphereBulkModulus	bulk modulus of sphere material			(common Gpa just be consistent)
	sphereShearModulus      modulus of sphere material			(common Gpa just be consistent)
	effectivePressure	confining pressure - n(pore pressure)		(to be consistent Gpa)
	alpha                   exponent in HM theory normally 1/3              (normally 1/3 suggested to be up to 1/6)
	"""
	C = aveNumContactPerSphere
	Phi = porosity
	K = sphereBulkModulus
	Mu = sphereShearModulus
	Pe = effectivePressure
	Nu = getPoissonRatioFromBulkAndShearModulus( K,Mu )		# compute poisson's ratio
	
	A0 = C * (1 - Phi) * Mu
	A = A0*A0 * Pe
	B0 = pi * (1 - Nu)
	B = 18.0 * B0*B0
	
	Khm = (A/B)**alpha
	
	return Khm

def hertzMindlinSpherePackEffectiveShearModulus( interGrainFriction,
						 aveNumContactPerSphere,
						 porosity, 
						 sphereBulkModulus, 
						 sphereShearModulus, 
						 effectivePressure,
						 alpha = 1.0/3.0 ):
	"""
	The function calculated the effective bulk modulus (Khm) of a random packing of spheres
	with average number of contacts per sphere (C), inter-sphere porosity (Phi), the shear
	and bulk of the spheres (Mu and K) along with the effective pressure (P),
	
	(Input Variable)	(Brief description)				(Units/Values)
	interGrainFriction      describes the adhesion between grains           (0=no friction between grains, 1=perfect adhesio )
	aveNumContactPerSphere	average number of contacts per sphere		(between 5 and 9 is common)
	porosity		inter-sphere space as fraction of whole		(fraction V/V normally around .36 dense random pack)
	sphereBulkModulus	bulk modulus of sphere material			(common Gpa just be consistent)
	sphereShearModulus      modulus of sphere material			(common Gpa just be consistent)
	effectivePressure	confining pressure - n(pore pressure)		(to be consistent Gpa)
	alpha                   exponent in HM theory normally 1/3              (normally 1/3 suggested to be up to 1/6)

	(Output Variable)	(Brief description)				(Units/Values)
	Mu_hm                   HM effective shear modulus of sphere pack       (consistent with input Gpa)
	"""
	f = interGrainFriction
	C = aveNumContactPerSphere
	Phi = porosity
	K = sphereBulkModulus
	Mu = sphereShearModulus
	Pe = effectivePressure
	Nu = getPoissonRatioFromBulkAndShearModulus( K,Mu )		# compute poisson's ratio

	A0 = 2.0 + 3.0*f - Nu*(1.0 + 3.0*f)
	A1 = 5.0*(2.0 - Nu)
	A = A0/A1
	B0 = C*(1.0 - Phi)*Mu
	B1 = 3.0 * B0*B0 * Pe
	B2 = pi*(1 - Nu)
	B = (B1/(2.0 * B2*B2))**alpha
	Mu_hm = A*B

	return Mu_hm

def dteModifiedHertzMindlinSpherePackEffectiveShearModulus( aveNumContactPerSphere,
							    porosity, 
							    sphereBulkModulus, 
							    sphereShearModulus, 
							    effectivePressure,
							    minEffectivePressure,
							    maxEffectivePressure,
							    maxBulkModulus,
							    accelerator = 1.0,
							    alpha = 1.0/3.0 ):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implementation allow the user to specify a
	maximum (dry-frame) bulk modulus associated with a min and max effective
	pressure range. The bulk modulus range is acheived by adjustment of the standard
	parameters to get the lower limit (by the user) then the upper limit is acheived by
	adjusment of the exponent (normally 1/3) to a value that reaches the upper bulk
	modulus limit. You can specify a modifier (accelerator) to more linearly sample the Khm 
	values, one suggestion would be to set it to the inverse ofupper limit of the power
	we compute. It maybe advantagious to compute Amin to a user specified Kmin at
	Pmin at a later date, currently I want the lower limit to be set by the value
	expected from the Hertz Mindlin theory

	
	(Input Variable)	(Brief description)				(Units/Values)
	aveNumContactPerSphere	average number of contacts per sphere		(between 5 and 9 is common)
	porosity		inter-sphere space as fraction of whole		(fraction V/V normally around .36 dense random pack)
	sphereBulkModulus	bulk modulus of sphere material			(common Gpa just be consistent)
	sphereShearModulus      modulus of sphere material			(common Gpa just be consistent)
	effectivePressure	confining pressure - n(pore pressure)		(to be consistent Gpa)
	minEffectivePressure    minimum effective presure to model              (Gpa)
	maxEffectivePressure    maximum effective presure to model              (Gpa)
	maxBulkModulus          maximum K for maximum pressure desired          (Gpa)
	accelerator             this control how fast the 
	alpha                   exponent in HM theory normally 1/3              (normally 1/3 probabaly leave this alone)


	(Output Variable)	(Brief description)				(Units/Values)
	Khm                     HM effective shear modulus of sphere pack       (consistent with input Gpa)
	"""
	
	C = aveNumContactPerSphere
	Phi = porosity
	K = sphereBulkModulus
	Mu = sphereShearModulus
	Pe = effectivePressure
	Nu = getPoissonRatioFromBulkAndShearModulus( K,Mu )		# compute poisson's ratio
	Pmax = maxEffectivePressure
	Pmin = minEffectivePressure
	Kmax = maxBulkModulus
	Ac   = accelerator

	
	A0 = (C * (1 - Phi) * Mu)/(sqrt(18.0)*pi * (1-Nu))
	A = A0*A0

	Amin = alpha
	Amax = log(Kmax)/log(A * Pmax)

	g = (Pe - Pmin) / (Pmax - Pmin)

	Anew = Amax * g**Ac + Amin * (1-g)**Ac
	
	Khm = (A * Pe)**Anew
	
	return Khm

def dteSecondModifiedHertzMindlinSpherePackEffectiveShearModulus( aveNumContactPerSphere,
								  porosity, 
								  sphereBulkModulus, 
								  sphereShearModulus, 
								  effectivePressure,
								  minEffectivePressure,
								  maxEffectivePressure,
								  maxBulkModulus,
								  accelerator = 1.0,
								  alpha = 1.0/3.0 ):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implementation allow the user to specify a
	maximum (dry-frame) bulk modulus associated with a min and max effective
	pressure range. The bulk modulus range is acheived by adjustment of the standard
	parameters to get the lower limit (by the user) then the upper limit is acheived by
	adjusment of the exponent (normally 1/3) to a value that reaches the upper bulk
	modulus limit. You can specify a modifier (accelerator) to more linearly sample the Khm 
	values, one suggestion would be to set it to the inverse ofupper limit of the power
	we compute. It maybe advantagious to compute Amin to a user specified Kmin at
	Pmin at a later date, currently I want the lower limit to be set by the value
	expected from the Hertz Mindlin theory

	
	(Input Variable)	(Brief description)				(Units/Values)
	aveNumContactPerSphere	average number of contacts per sphere		(between 5 and 9 is common)
	porosity		inter-sphere space as fraction of whole		(fraction V/V normally around .36 dense random pack)
	sphereBulkModulus	bulk modulus of sphere material			(common Gpa just be consistent)
	sphereShearModulus      modulus of sphere material			(common Gpa just be consistent)
	effectivePressure	(litho-static pressure) - n(pore pressure)	(to be consistent Gpa)
	minEffectivePressure    minimum effective presure to model              (Gpa)
	maxEffectivePressure    maximum effective presure to model              (Gpa)
	maxBulkModulus          maximum K for maximum pressure desired          (Gpa)
	accelerator             this control how fast the 
	alpha                   exponent in HM theory normally 1/3              (normally 1/3 probabaly leave this alone)


	(Output Variable)	(Brief description)				(Units/Values)
	Khm                     HM effective shear modulus of sphere pack       (consistent with input Gpa)
	"""
	
	C = aveNumContactPerSphere
	Phi = porosity
	K = sphereBulkModulus
	Mu = sphereShearModulus
	Pe = effectivePressure
	Nu = getPoissonRatioFromBulkAndShearModulus( K,Mu )		# compute poisson's ratio
	Pmax = maxEffectivePressure
	Pmin = minEffectivePressure
	Kmax = maxBulkModulus
	Ac   = accelerator

	
	A0 = (C * (1 - Phi) * Mu)/(sqrt(18.0)*pi * (1-Nu))
	A = A0*A0

	KcalcMax = (A * Pmax)**alpha
	e = (Kmax / KcalcMax)**(1.0/alpha)

	g = (Pe - Pmin) / (Pmax - Pmin)

	m = 1.0 + (e-1.0)*g
	
	Khm = (A * m*Pe)**alpha
	
	return Khm


def alphaHM( aveNumContactPerSphere,
	     porosity, 
	     sphereBulkModulus, 
	     sphereShearModulus, 
	     effectivePressure):
	"""
	This function computes the alpha constant in the Hertz-Mindlin equation if we
	parameterize the equation as Khm = alpha * P**(1/3) which means:

	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3)

	In my modifications this can be assumed to be a fitting parameter.

	"""

	C = aveNumContactPerSphere
	Phi = porosity
	K = sphereBulkModulus
	Mu = sphereShearModulus
	Pe = effectivePressure

	Nu = (3.0*K - 2.0*Mu)/(2.0*(3.0*K + Mu))

	a0 = (C * (1-Phi) * Mu) / (sqrt(18.0) * pi * (1-Nu))

	alpha = a0**(2.0/3.0)

	return alpha

def dteThirdModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								firstRefPressure,
								firstRefBulkModulus,
								secondRefPressure,
								secondRefBulkModulus):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has two reference
	points (P1,K1) and (P2,K2), these effective pressure and bulk mudulus pairs will
	be used in the following mammer to generate modified Hertz-Mindlin points. First
	I assume a useful modified form would be a shifted scaled version of the original
	Hertz-Mindlin equation. As usual, the assumed parametrization of the Hertz-Mindlin
	equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	I make the modification:

	     KhmMod = alpha' * (P + dP)***1/3)

	where

	     alpha' and dP are assumed to be variables to be fitted to our input data.

	The process first solves for dP given by:

	     dp = ((K1/K2)^3 * P2 - P1)/(1 - (K1/K2)^3)

	Then we solve for alpha' from either one of the input conditions

	     alpha' = K1/(P1 + dP)^(1/3) = K2/(P1 + dP)^(1/3)
	     
	This function then uses the modified equation to calculte KhmMod at a new
	pressure
	
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Gpa for consistency with moduli)
	firstRefPressure	Pressure for first reference K and alpha	(Gpa for consistency with moduli)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)
	secondRefPressure	Pressure for second reference K and alpha	(Gpa for consistency with moduli)
	secondRefBulkModulus    Second reference bulk modulus			(Gpa)


	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""
	P = newPressure
	P1 = firstRefPressure
	K1 = firstRefBulkModulus
	P2 = secondRefPressure
	K2 = secondRefBulkModulus

	K3 = (K1/K2)**3
	dp = (K3*P2 - P1) / (1.0 - K3)      # Compute the pressure shift term

	if (P + dp) <= 0.0:
		
		KhmMod = 0
		return KhmMod

	else:

		Ap = K1/((P1 + dp)**(1.0/3.0))      # Compute alpha'
		KhmMod = Ap * (P + dp)**(1.0/3.0)
	
	return KhmMod

def dteThirdModifiedHertzMindlinSpherePackEffectiveShearModulus( KhmMod, beta=-999.25, nu = 0.094, f = 1.0 ):
	"""
	This function assumes that you already have a list/set of modified Khm (perhaps from a run of
	'dteThirdModifiedHertzMindlinSpherePackEffectiveBulkModulus' which is then going to be modified
	to arrive at the modified Hertz-Mindlin shear modulus. This makes the observation that:

        	MuHm = beta*Khm
		
	where (if not specifically specified by user):
	
	        beta = 3*(2 + 3*f - nu*(1 + 3*f))) / (5*(2 - nu)
		     f = grain contact adhesion factor
		     nu = Poisson's ratio of the grains
		     
	other wise, it is the value input by user.

	Note: the input KhmMod is assumed to be a list, so make sure it is.
	The output  will be a list with the associated MuHmMod values
	
	"""
	MuHmMod = []
	
	if beta == -999.25:
		
		beta = 3.0 * (2.0 + 3.0*f - nu*(1.0 + 3.0*f)) / (5.0*(2.0 - nu))
		MuHmMod = [beta * k for k in KhmMod]
		
	else:
		
		MuHmMod = [beta * k for k in KhmMod]

	return MuHmMod

def findShiftFunctionMaker( k1, p1, k2, p2 ):
	"""
	This function is used in my fourth modification of the Hertz-Mindlin theory.
	The input to this function is just the reference bulk modulus and effective
	pressure pairs at the reference porosity to be used. This function then
	creates a function to be used to find the shift term for modification (IV)

	"""
	return lambda x: (k1-k2)*x**(1.0/3.0) + k2*(p1+x)**(1.0/3.0) - k1*(p2+x)**(1.0/3.0)


def dteFourthModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								 firstRefPressure,
								 firstRefBulkModulus,
								 secondRefPressure,
								 secondRefBulkModulus):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has two reference
	points (P1,K1) and (P2,K2), these effective pressure and bulk mudulus pairs will
	be used in the following mammer to generate modified Hertz-Mindlin points. First
	I assume a useful modified form would be a shifted scaled version of the original
	Hertz-Mindlin equation. As usual, the assumed parametrization of the Hertz-Mindlin
	equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	I make the modification:

	     KhmMod = alpha' * [(P + dP)***1/3 - dP**1/3]

	where

	     alpha' and dP are assumed to be variables to be fitted to our input data.

	This modification insures KhmMod = 0 when P = 0. Again we will solve for dP
	within the non-linear equation below numerically

	     f(dP) = (K1-K2)*dP**(1/3) + K2*(P1+dP)**(1/3) - K1*(P2+dP)**(1/3).

	In other words, we are finding the root of the equation above numerically f(dP) = 0.

	Then we solve for alpha' from either one of the input conditions

	     alpha' = K1/((P1 + dP)^(1/3) - dP**(1/3)) = K2/((P1 + dP)^(1/3) - dP**(1/3))
	     
	This function then uses the modified equation to calculte KhmMod at a new
	pressure
	
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Gpa for consistency with moduli)
	firstRefPressure	Pressure for first reference K and alpha	(Gpa for consistency with moduli)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)
	secondRefPressure	Pressure for second reference K and alpha	(Gpa for consistency with moduli)
	secondRefBulkModulus    Second reference bulk modulus			(Gpa)


	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""
	P = newPressure
	P1 = firstRefPressure
	K1 = firstRefBulkModulus
	P2 = secondRefPressure
	K2 = secondRefBulkModulus

	f = findShiftFunctionMaker( K1, P1, K2, P2 )           # Create function from which to find roots
	dp = float(fsolve(f,1.0)[0])                                     # find root starting from 1.0
#	print "dp = ",dp

	Ap = K1/((P1 + dp)**(1.0/3.0) - (dp)**(1.0/3.0))       # Compute alpha'
#	print "alpha = ",float(Ap)
	KhmMod = Ap * ((P + dp)**(1.0/3.0) - (dp)**(1.0/3.0))  # Compute modified Khm
	
	return KhmMod

def dteFifthModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								RefPressure,
								RefBulkModulus,
								scale = 1.0,
								beta = (1.0/3.0)):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has two reference
	points (P1,K1) and (P2,K2), these effective pressure and bulk mudulus pairs will
	be used in the following mammer to generate modified Hertz-Mindlin points. First
	I assume a useful modified form would be a shifted scaled version of the original
	Hertz-Mindlin equation. As usual, the assumed parametrization of the Hertz-Mindlin
	equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	I make the modification:

	     KhmMod = alpha' * P**1/3

	where

	     alpha' is assumed to be variable and fitted to our input data.
	     alpha' = Kref/Pref**(1/3)
	     
	This function then uses the modified equation to calculte KhmMod at a new
	pressure
	
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Gpa for consistency with moduli)
	firstRefPressure	Pressure for first reference K and alpha	(Gpa for consistency with moduli)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)
	secondRefPressure	Pressure for second reference K and alpha	(Gpa for consistency with moduli)
	secondRefBulkModulus    Second reference bulk modulus			(Gpa)


	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""
	P  = newPressure
	P1 = RefPressure
	K1 = RefBulkModulus

	Ap = K1 / (P1**beta)        # Compute alpha'

	KhmMod = K1 * (P/P1)**(scale*beta)       # Compute modified Khm
	
	return KhmMod

def dteFifthModifiedHertzMindlinSpherePackEffectiveShearModulus( KhmMod, MuRef, Kref, delta=-999.25):
	"""
	This function assumes that you already have a list/set of modified Khm (from a run of
	'dteFifthModifiedHertzMindlinSpherePackEffectiveBulkModulus' which is then going to be modified
	to arrive at the modified Hertz-Mindlin shear modulus. This makes the observation that:

        	MuHm = beta*Khm

	however, we want to have MuHm = MuRef if Khm = Kref. This can be acheived if we modify the
	equation above as:

	        MuHmMod = MuRef * (Khm/KRef)**delta

	Where delta plays the role of a stretch factor if not specified is set to one, otherwise,
	it will be set to the user value.

	Note: the input KhmMod is assumed to be a list, so make sure it is.
	The output  will be a list with the associated MuHmMod values
	
	"""
	MuHmMod = []
	
	if delta == -999.25:
		
		delta = 1.0
		MuHmMod = [MuRef * (k/Kref)**delta for k in KhmMod]
		
	else:
		
		MuHmMod = [MuRef * (k/Kref)**delta for k in KhmMod]

	return MuHmMod

def dteSixthModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								RefPressure,
								RefBulkModulus,
								scale = 1.0,
								pow = (1.0/3.0)):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has two reference
	points (P1,K1) and (P2,K2), these effective pressure and bulk mudulus pairs will
	be used in the following mammer to generate modified Hertz-Mindlin points. First
	I assume a useful modified form would be scaled version of the original, with the
	reference point given by its effective pressure and bulk modulus uneffected. As
	usual, the assumed parametrization of the Hertz-Mindlin equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	I make the modification:

	     KhmMod =  scale * alpha * P**1/3 + (1-scale) * Kref

	where

	     alpha = Kref/Pref**(1/3)

	This insures the KhmMod = Kref when P = Pref and the scale factor behaves a:

	     scale > 1 => stretching
	     scale = 1 => no change
	     scale < 1 => squeezing

	This function then uses the modified equation to calculte KhmMod at a new
	pressure. It should be noted some of your K's may become negative from the
	stretching which can either be set to zero or used as is since it will be the
	deltas we are interested in and in the region of the real data it is probably
	well behaved.
	
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Mpa)
	firstRefPressure	Pressure for first reference K and alpha	(Mpa)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)

	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""
	P  = newPressure * .001                                   # Conversion from Mpa to Gpa
	P1 = RefPressure * .001                                   # Conversion from Mpa to Gpa
	K1 = RefBulkModulus

	Ap = K1 / (P1**pow)                                       # Compute alpha

	KhmMod = scale * Ap * (P**pow) + (1.0-scale) * K1         # Compute modified Khm
	
	return KhmMod

def dteSixthModifiedHertzMindlinSpherePackEffectiveShearModulus( KhmMod, MuRef, Kref, scale = 1.0 ):
	"""
	This function assumes that you already have a list/set of modified Khm (from a run of
	'dteSixthModifiedHertzMindlinSpherePackEffectiveBulkModulus' which is then going to be modified
	to arrive at the modified Hertz-Mindlin shear modulus. This makes the observation that:

        	MuHm = beta*Khm

	however, we want to have MuHm = MuRef if Khm = Kref. This can be acheived if we modify the
	equation above as:

	        MuHmMod = (MuRef/Kref) * Khm

	The scaling work much as in the bulk modulus version on this modified equation making the
	reference point stationary throughout the stretching/sqeezing:

	        MuHmMod = scale * (MuRef/Kref) * Khm + (1-scale) * MuRef

	Note: the input KhmMod is assumed to be a list, so make sure it is.
	The output  will be a list with the associated MuHmMod values
	
	"""
	MuHmMod = []
	
	MuHmMod = [scale * (MuRef/Kref) * k**pow + (1.0-scale) * MuRef for k in KhmMod]

	return MuHmMod

def dteSeventhModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								  RefPressure,
								  RefBulkModulus,
								  scale = 1.0,
								  expo = (1.0/3.0)):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has a reference
	point (P0,K0) this effective pressure and bulk mudulus pair will be used in
	the following mammer to generate modified Hertz-Mindlin points. First, I
	assume a useful modified form would be a scaled version of the original, with
	the reference point given by its effective pressure and bulk modulus uneffected.
	As usual, the assumed parametrization of the Hertz-Mindlin equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	In this implimentation we will assume that alpha is completely defined by the
	reference point (P0,K0). We make the modification:

	     KhmMod =  alpha * (scale * P)**1/3 - (alpha * (scale * P0)**1/3)  + K0
	            =  alpha * scale**1/3 * (P**1/3 - P0**1/3) +  K0

	where

	     alpha = K0/P0**(1/3)

	This insures the KhmMod = Kref when P = Pref and the scale factor behaves a:

	     scale > 1 => stretching
	     scale = 1 => no change
	     scale < 1 => squeezing

	This function then uses the modified equation to calculte KhmMod at a new
	pressure. It should be noted some of your K's may become negative from the
	stretching which can either be set to zero or used as is since it will be the
	deltas we are interested in for 4D type analysis and in the region of real
	data it is probably well behaved.
	
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Mpa)
	firstRefPressure	Pressure for first reference K and alpha	(Mpa)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)

	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""
	P  = float(newPressure) * .001                                   # Conversion from Mpa to Gpa
	P0 = float(RefPressure) * .001                                   # Conversion from Mpa to Gpa
	K0 = float(RefBulkModulus)
	s = float(scale)**float(expo)

	alpha = K0 / (P0**expo)                                      # Compute alpha

	KhmMod = alpha * s * (P**expo - P0**expo) +  K0         # Compute modified Khm
	
	return KhmMod

def dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus( KhmMod, MuRef, Kref, scaleMu = 1.0, scaleK = 1.0, expo = 1.0/3.0, **kwargs ):
	"""
	This function assumes that you already have a list/value of modified Khm (from a run of
	'dteSeventhModifiedHertzMindlinSpherePackEffectiveBulkModulus' which is then going to be modified
	to arrive at the modified Hertz-Mindlin shear modulus. This makes the observation that:

        	MuHm = beta*Khm


	however, we want to have MuHm = MuRef if Khm = Kref. This can be acheived if we modify the
	equation above as:

	        MuHmMod = (MuRef/Kref) * Khm

	The scaling works much as in the bulk modulus version on this modified equation making the
	reference point stationary throughout the stretching/sqeezing:

	        s = (scaleMu/scaleK)**1/3
		beta = (MuRef/Kref)
		MuHmMod = s * beta * Khm + (1.0-s)*Kref

	This in essence removes the effect of the K scalar (scaleK)  and replaces it with the Mu scalar scaleMu)

	Note: the input KhmMod is assumed to be a list, so make sure it is.
	The output  will be a list with the associated MuHmMod values
	
	"""
	sK = float(scaleK)
	sMu = float(scaleMu)
	s = (sMu/sK)**expo
	if isinstance(MuRef,(int, long, float, complex)) and isinstance(Kref,(int, long, float, complex)):
		beta = MuRef/Kref
	else:
		print " "
		print "**** Error: dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus ****"
		print "     Reference bulk and shear moduli (MuRef,Kref) must be scalars."
		print " "
		return
	
	if type( KhmMod ) == list:
		MuHmMod = []
		MuHmMod = [s * beta * k + (1.0 - s)*Kref  for k in KhmMod]
		
	elif isinstance( KhmMod,(int, long, float, complex)):
		MuHmMod = s * beta * KhmMod + (1.0 - s)*Kref

	elif type( KhmMod ) == np.ndarray:
		MuHmMod = s * beta * KhmMod + (1.0 -s)*Kref

	else:
		print " "
		print "**** Error: dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus ****"
		print "     Input bulk modulus is not of a datatype recognized by this function."
		print " "
		return

	return MuHmMod

def dteEighthModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								 RefPressure,
								 RefBulkModulus,
								 scale = 1.0,
								 blendPow = 3.0, 
								 expo = (1.0/3.0)):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has a reference
	point (P0,K0) this effective pressure and bulk mudulus pair will be used in
	the following mammer to generate modified Hertz-Mindlin points. First, I
	assume a useful modified form would be a scaled version of the original, with
	the reference point given by its effective pressure and bulk modulus uneffected.
	As usual, the assumed parametrization of the Hertz-Mindlin equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	In this implimentation we will assume that alpha is completely defined by the
	reference point (P0,K0). We make the modification:

	     KhmMod =  alpha * (scale * P)**1/3 - (alpha * (scale * P0)**1/3)  + K0
	            =  alpha * scale**1/3 * (P**1/3 - P0**1/3) +  K0

	where

	     alpha = K0/P0**(1/3)

	This insures the KhmMod = Kref when P = Pref and the scale factor behaves as follows:

	     scale > 1 => stretching
	     scale = 1 => no change
	     scale < 1 => squeezing

	There is a problem with this implementation at lower pressures where the computed
	K's may become negative. To compensate for this effect, we can set the near zero
	effective pressure behaviour to be like the original formulation, where:

	     Khm1 = alpha * P**(1/3)

	and for effective pressure equal to or greater than P0 the behave as:

	     Khm2 =  alpha * scale**1/3 * (P**1/3 - P0**1/3) +  K0.

	If you examine Khm2 it will become apparent that as scale -> 1 => Khm2 -> Khm1
	Using this fact and letting:

	     s = scale**1/3, e = (P/P0)**1/3 and sub back alpha = K0/P0**(1/3)
	     
	we rewrite Khm2 as:

	     Khm2 = [s*e + (1-s)]*K0, again note Khm2 -> Khm1 when s -> 1.

	now for 0 < P < P0 we use this modification:

	     Khm2' = [s'*e + (1-s')]*K0

	where:

	     s' = 1 + (s - 1) * e**d.

	As can be seen if e=1 s'=s and if e=0 s'=1 and e=1 when P = P0 and e=0 when P = 0.
	This is exactly what we want. "d" is a power we want to choose so that Khm2' > 0 for
	0<P<P0
	     

	
	This function then uses the modified equation to calculate KhmMod at a new
	pressure. 
	
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Mpa)
	firstRefPressure	Pressure for first reference K and alpha	(Mpa)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)
	scale                   Scaling parameter to match data                 (None)
	blendPow                blending variable to determine                  (None)
	                        behavior at lower pressures
	expo                    Hertz Mindlin exponent (standard 1/3)           (None)
	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""

	d = blendPow
	P  = float(newPressure) * .001                                   # Conversion from Mpa to Gpa
	P0 = float(RefPressure) * .001                                   # Conversion from Mpa to Gpa
	K0 = float(RefBulkModulus)
	s = float(scale)**float(expo)

	e = (P/P0)**float(expo)

	sp = 1.0 + (s - 1.0) * e**d


	
	if P >= 0.0 and P < P0:
		KhmMod =  ( sp*e + (1.0 -sp) ) * K0         # Alternate fomulation for 0 < P < 1
		
	elif P >= P0:
		KhmMod = ( s*e + (1.0 - s) ) * K0           # Same as 7th in different form

	else:
		print " "
		print "**** Error: dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus ****"
		print "            Effective pressure input cannot be negative."
		print " "
		return
	
	return KhmMod

def dteEighthModifiedHertzMindlinSpherePackEffectiveShearModulus( KhmMod, MuRef, Kref, **kwargs ):
	"""
	(Note: this is identical to the seventh modification except the scaling will not be changed)
	This function assumes that you already have a list/value of modified Khm (from a run of
	'dteSeventhModifiedHertzMindlinSpherePackEffectiveBulkModulus' which is then going to be modified
	to arrive at the modified Hertz-Mindlin shear modulus. This makes the observation that:

        	MuHm = beta*Khm


	however, we want to have MuHm = MuRef if Khm = Kref. This can be acheived if we modify the
	equation above as:

	        MuHmMod = (MuRef/Kref) * Khm

	where:
		beta = (MuRef/Kref)

	"""

	if isinstance(MuRef,(int, long, float, complex)) and isinstance(Kref,(int, long, float, complex)):
		beta = MuRef/Kref
	else:
		print " "
		print "**** Error: dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus ****"
		print "     Reference bulk and shear moduli (MuRef,Kref) must be scalars."
		print " "
		return
	
	if type( KhmMod ) == list:
		MuHmMod = [beta * k  for k in KhmMod]
		
	elif isinstance( KhmMod,(int, long, float, complex)):
		MuHmMod = beta * KhmMod

	elif type( KhmMod ) == np.ndarray:
		MuHmMod = beta * KhmMod

	else:
		print " "
		print "**** Error: dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus ****"
		print "     Input bulk modulus is not of a datatype recognized by this function."
		print " "
		return

	return MuHmMod

def dteNinthModifiedHertzMindlinSpherePackEffectiveBulkModulus( newPressure,
								RefPressure,
								RefBulkModulus,
								scale = 1.0,
								expo = (1.0/3.0)):
	"""
	This function is a modified version of the Hertz Mindlin random sphere pack
	bulk modulus function. This implimentation assumes the user has a reference
	point (P0,K0) this effective pressure and bulk mudulus pair will be used in
	the following mammer to generate modified Hertz-Mindlin points. First, I
	assume a useful modified form would be a scaled version of the original, with
	the reference point given by its effective pressure and bulk modulus uneffected.
	As usual, the assumed parametrization of the Hertz-Mindlin equation is:

	     Khm = alpha * P**(1/3)
	     
	where
	
	     alpha = ((C*(1-Phi)*Mu)/(sqrt(18)*pi*(1-Nu)))**(2/3).

	In this implimentation we will assume that alpha is completely defined by the
	reference point (P0,K0). We make the modification:

	     KhmMod =  K0 * (P/P0)**(scale/3)
	
	where

	     alpha = K0/P0**(1/3)

	This insures the KhmMod = Kref when P = Pref and the scale factor behaves as follows:

	     scale > 1 => stretching of pressure axis
	     scale = 1 => no change (simple Hertz-Mindlin formulation with fixed point (K0.P0))
	     scale < 1 => squeezing of pressure axis

	Another way to look at this is KhmMod = Khm when scale = 1. When scale > 1 the (KhmMod,P)
	plot is steeper compared to the standard Hertz-Mindlin formulation. When scale < 1 the
	(KhmMod,P) plot is less steep compared to the standard Hertz-Mindlin formulation.
	
	(Input Variable)	(Brief description)				(Units/Values)
	newPressure           	Pressure to compute new modified HM K		(Mpa)
	firstRefPressure	Pressure for first reference K and alpha	(Mpa)
	firstRefBulkModulus     First reference bulk modulus			(Gpa)
	scale                   Scaling parameter to match data                 (None)

	(Output Variable)	(Brief description)				(Units/Values)
	KhmMod                  Modified HM bulk modulus                        (consistent with input Gpa)
	"""

	P  = float(newPressure) * .001                                   # Conversion from Mpa to Gpa
	P0 = float(RefPressure) * .001                                   # Conversion from Mpa to Gpa
	K0 = float(RefBulkModulus)
	s = float(scale)
	
	if P < 0.0:
		print " "
		print "**** Error: dteSeventhModifiedHertzMindlinSpherePackEffectiveShearModulus ****"
		print "            Effective pressure input cannot be negative."
		print " "
		return
	
	e = (P/P0)**float(expo)

	KhmMod = K0 * e**s

	return KhmMod

def betaFunction( K0, Mu0, Kmatrix, MuMatrix, Phi0, Phi1, beta=1.472 ):
	"""
	This function is used by "findBeta" to iterate for a beta value in the
	equation: Mu = beta * K so that the soft sand Mu HS lower bound will
	give Mu0 at Phi0

	"""
	
	K1p = softSandEffectiveBulkModulusHertzMindlinAtNewPorosity( Phi0, Phi1, K0, Mu0, Kmatrix )
	Mu1p = beta * K1p
	Mu0p = softSandEffectiveShearModulusHertzMindlinAtNewPorosity( Phi1, Phi0, K1p, Mu1p, MuMatrix )
	delta = Mu0 - Mu0p

	return delta

def findBeta( K0, Mu0, Kmatrix, MuMatrix, Phi0, Phi1 ):
	"""
	This function finds 'beta' that is used to compute MuHm = beta * Khm. The beta is
	ment to insure that when we use the soft sand HS lower bound for Mu we get Mu0
	for Phi0. The root is found using the bisection method implimented in scipy.

	"""

	betaF = lambda x : betaFunction(K0, Mu0, Kmatrix, MuMatrix, Phi0, Phi1, x)
	beta = bisect( betaF, 0.2, 3.0 )

	return beta

def findHertzMindlinModuliFromReferenceState( Phi1, K0, Mu0, Phi0, Kmatrix, MuMatrix, **kwargs):
	"""
	This function is used by 'dryFrameModuliChangeWithPressure' to find the Hertz-Mindlin
	reference (KhmRef,MuHmRef) at Phi1 from the reference point (K0,Mu0) at Phi0. This is
	done at constant effective pressure P0. The units are consistent with the calling
	function 'dryFrameModuliChangeWithPressure', so look there for description.

	"""

	# define a two equation functionfor numerically finding two roots associated with input variables 
	f = lambda x : [softSandEffectiveBulkModulusHertzMindlinAtNewPorosity(Phi1, Phi0, x[0], x[1], Kmatrix) - K0,
			softSandEffectiveShearModulusHertzMindlinAtNewPorosity(Phi1, Phi0, x[0], x[1], MuMatrix) - Mu0]
	
	# establish starting point for root search
	Ki = softSandEffectiveBulkModulusHertzMindlinAtNewPorosity(Phi0, Phi1, K0, Mu0, Kmatrix)
	MuI = softSandEffectiveShearModulusHertzMindlinAtNewPorosity(Phi0, Phi1, K0, Mu0, MuMatrix)

	# Use scipy routine 'fsolve' to find the roots
	[KhmP0, MuHmP0] = list( fsolve(f,[Ki,MuI]) )

	return KhmP0, MuHmP0
	

def dryFrameModuliChangeWithPressure( Pnew, Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, sigma = 1.0, scaleK=1.0, blendPow=1.0, scaleMu=1.0 , printVar=0, **kwargs):
	"""
	This funcition uses the modified Hertz-Mindlin method to update the original
	dry-frame bulk modulus under a change in effective pressure. Matrix properties
	(Kmatrix, MuMatrix) are assumed to be unchanged. The reference state or original
	state is given by (K0,Mu0,Phi0,P0) = reference(bulk modulus, shear modulus,
	porosity, effective pressure) the output will be (K1,Mu1) = resultant(bulk modulus,
	shear modulus) at the new effective pressure P1. The assumption that an intermediate
	state, favoured by the HM theory, is assumed to occure at porosity of 0.36 (V/V). Though
	the predicted HS K and Mu's can be independently stretch squeezed, in practice, it
	should only be the HS K's that is squeezed and use the resulting HS Mu's as is; however,
	it maybe necessary to scale the Mu's independently to fit the actual data you may have.

	(Input Variables)

	(Variable)    (Description)                                 (Unit)
	Pnew          New pressure to compute K1,Mu1                (Mpa)
	Pref          effective pressure at reference state         (Mpa)
	Kref          dryframe bulk modulus at reference state      (Gpa)
	MuRef         shear modulus at reference state              (Gpa)
	PhiRef        porosity at reference state                   (V/V)
	Kmatrix       Matrix minerals net bulk modulus              (Gpa)
	MuMatrix      Matrix minerals net bulk modulus              (Gpa)
	sigma         scaling factor for 9th mod HM
	scaleK        scaling factor to stretch/squeeze Khm  (Not used in ninth version of HM
	scaleMu       scaling factor to stretch/squeeze MuHm (Not used in eighth version of HM
	              scaleK = scaleMu)
	              resulting K's and Mu's for pressure range
		                 If scale > 1 => stretching
		                 If scale = 1 no stretch/squeeze
				 If scale < 1 => squeezing

	(Output Variables)
	K1            dry-frame bulk modulus at Pnew                (Gpa)
	Mu1           dry-frame shear modulus at Pnew               (Gpa)
	"""
	s = float(sigma)
	sK = float(scaleK)
	sMu = float(scaleMu)
	bp = float(blendPow)

	Phi0 = float(PhiRef)
	Phi1 = 0.36

	P1 = float(Pnew)*0.001           # Conversion from Mpa to Gpa
	P0 = float(Pref)*0.001

	K0 = float(Kref)
	Mu0 = float(MuRef)

	Kmat = float(Kmatrix)
	MuMat = float(MuMatrix)
	
	KhmP0, MuHmP0 = findHertzMindlinModuliFromReferenceState(Phi1, K0, Mu0, Phi0, Kmat, MuMat)


	KhmP1 = dteNinthModifiedHertzMindlinSpherePackEffectiveBulkModulus( P1, P0, KhmP0, s )
	MuHmP1 = dteEighthModifiedHertzMindlinSpherePackEffectiveShearModulus( KhmP1, Mu0, K0 )

	if printVar == 1:
			print " "
			print "Hertz Mindlin reference K and Mu at Phi=0.36 and P="+str(P0*1000.0)+" :"+str(KhmP0)+" , "+str(MuHmP0)
			print "Hertz Mindlin K and Mu at Phi="+str(Phi1)+" and P="+str(P1*1000.0)+" :"+str(KhmP1)+" , "+str(MuHmP1)
			print " "

	K1 = softSandEffectiveBulkModulusHertzMindlinAtNewPorosity(Phi1, Phi0, KhmP1, MuHmP1, Kmatrix)
	Mu1 = softSandEffectiveShearModulusHertzMindlinAtNewPorosity(Phi1, Phi0, KhmP1, MuHmP1, MuMatrix)
	

	return K1, Mu1
	
	
	
	

def softSandEffectiveBulkModulusHertzMindlinAtNewPorosity( originalPorosity,
							   newPorosity,
							   bulkModulusHM,
							   shearModulusHM,
							   grainBulkModulus ):
	"""
	Using a heuristically modified Hashin-Shtrikman lower bound to adjust the Hertz-Mindlin
	derived bulk modulus (Khm) and shear modulus (MuHm) to effective moduli (Keff and MuEff)
	at a different porosity.

	(Input Variable)	(Brief description)				(Units/Values)
	originalPorosity        porosity used in Hertz-Mindlin calculation      (fraction V/V)
	newPorosity             porosity for effective porosity                 (fraction V/V)
	bulkModulusHM           Hertz-Mindlin bulk modulus                      (common Gpa just be consistent)
	shearModulusHM          Hertz-Mindlin shear modulus                     (common Gpa just be consistent)
	grainBulkModulus        sphere bulk modulus in Hertz-Mindlin calc       (common Gpa just be consistent)
	
	(Output Variable)	(Brief description)				(Units/Values)
	K_eff_soft              new effective bulk modulus                      (common Gpa just be consistent)
	"""

	Phi0 = originalPorosity
	Phi = newPorosity
	Khm = bulkModulusHM
	Mu_hm = shearModulusHM
	Ks = grainBulkModulus
	
	PhiR = Phi/Phi0
	
	if Khm == 0.0:
		K_eff_soft = 0.0
	else:
		A0 = PhiR/(Khm + 4.0*Mu_hm/3.0)
		A1 = (1.0-PhiR)/(Ks + 4.0*Mu_hm/3.0)
		A = 1.0 / (A0 + A1)
		B = 4.0*Mu_hm/3.0
		
		K_eff_soft = A - B

	return K_eff_soft


def softSandEffectiveShearModulusHertzMindlinAtNewPorosity( originalPorosity,
							    newPorosity,
							    bulkModulusHM,
							    shearModulusHM,
							    grainShearModulus ):
								       
	"""
	Using a heuristically modified Hashin-Shtrikman lower bound to adjust the Hertz-Mindlin
	derived bulk modulus (Khm) and shear modulus (MuHm) to effective moduli (Keff and MuEff)
	at a different porosity.

	(Input Variable)	(Brief description)				(Units/Values)
	originalPorosity        porosity used in Hertz-Mindlin calculation      (fraction V/V)
	newPorosity             porosity for effective porosity                 (fraction V/V)
	bulkModulusHM           Hertz-Mindlin bulk modulus                      (common Gpa just be consistent)
	shearModulusHM          Hertz-Mindlin shear modulus                     (common Gpa just be consistent)
	grainShearModulus       sphere shear modulus inHertz-Mindlin calc       (common Gpa just be consistent)
	
	(Output Variable)	(Brief description)				(Units/Values)
	Mu_eff_soft             new effective shear modulus                     (common Gpa just be consistent)
	"""
	Phi0 = originalPorosity
	Phi = newPorosity
	Khm = bulkModulusHM
	Mu_hm = shearModulusHM
	Mu_s = grainShearModulus
	
	PhiR = Phi/Phi0

	if Khm == 0.0:
		Mu_eff_soft = 0.0
	else:
		A = Mu_hm/6.0 * ((9.0*Khm + 8.0*Mu_hm) / (Khm + 2.0*Mu_hm))
		B0 = PhiR / (Mu_hm + A)
		B1 = (1.0 - PhiR) / (Mu_s + A)
		B = 1.0 / (B0 + B1)
		
		Mu_eff_soft = B - A

	return Mu_eff_soft

def stiffSandEffectiveBulkModulusHertzMindlinAtNewPorosity( originalPorosity,
							    newPorosity,
							    bulkModulusHM,
							    grainBulkModulus, 
							    grainShearModulus ):
	"""
	Using a heuristically modified Hashin-Shtrikman lower bound to adjust the Hertz-Mindlin
	derived bulk modulus (Khm) and shear modulus (MuHm) to effective moduli (Keff and MuEff)
	at a different porosity.

	(Input Variable)	(Brief description)				(Units/Values)
	originalPorosity        porosity used in Hertz-Mindlin calculation      (fraction V/V)MuHmP1
	newPorosity             porosity for effective porosity                 (fraction V/V)
	bulkModulusHM           Hertz-Mindlin bulk modulus                      (common Gpa just be consistent)
	grainBulkModulus        sphere bulk modulus in Hertz-Mindlin calc       (common Gpa just be consistent)
	grainShearModulus       sphere shear modulus in Hertz-Mindlin calc      (common Gpa just be consistent)
	
	(Output Variable)	(Brief description)				(Units/Values)
	K_eff_stiff             new effective bulk modulus                      (common Gpa just be consistent)
	"""
	Phi0 = originalPorosity
	Phi = newPorosity
	K_hm = bulkModulusHM
	K_s = grainBulkModulus
	Mu_s = grainShearModulus
	
	PhiR = Phi/Phi0
	
	if K_hm == 0.0:
		K_eff_stiff = 0.0
	else:
		A0 = PhiR / (K_hm + 4.0*Mu_s/3.0)
		A1 = (1-PhiR) / (K_s + 4.0*Mu_s/3.0)
		A = 1.0 / (A0 + A1)
		B = 4.0 * Mu_s / 3.0
		
		K_eff_stiff = A - B

	return K_eff_stiff


def stiffSandEffectiveShearModulusHertzMindlinAtNewPorosity( originalPorosity,
							     newPorosity,
							     shearModulusHM,
							     grainBulkModulus,
							     grainShearModulus ):
								       
	"""
	Using a heuristically modified Hashin-Shtrikman lower bound to adjust the Hertz-Mindlin
	derived bulk modulus (Khm) and shear modulus (MuHm) to effective moduli (Keff and MuEff)
	at a different porosity.

	(Input Variable)	(Brief description)				(Units/Values)
	originalPorosity        porosity used in Hertz-Mindlin calculation      (fraction V/V)
	newPorosity             porosity for effective porosity                 (fraction V/V)
	shearModulusHM          Hertz-Mindlin shear modulus                     (common Gpa just be consistent)
	grainBulkModulus        sphere shear modulus inHertz-Mindlin calc       (common Gpa just be consistent)
	grainShearModulus       sphere shear modulus inHertz-Mindlin calc       (common Gpa just be consistent)
	
	(Output Variable)	(Brief description)				(Units/Values)
	Mu_eff_stiff            new effective shear modulus                     (common Gpa just be consistent)
	"""
	Phi0 = originalPorosity
	Phi = newPorosity
	Mu_hm = shearModulusHM
	K_s = grainBulkModulus
	Mu_s = grainShearModulus
	
	PhiR = Phi/Phi0
	A = Mu_s/6.0 * ((9.0*K_s + 8.0*Mu_s) / (K_s + 2.0*Mu_s))
	B0 = PhiR / (Mu_hm + A)
	B1 = (1 - PhiR) / (Mu_s + A)
	B = 1.0 / (B0 + B1)
	
	Mu_eff_stiff = B - A

	return Mu_eff_stiff

def voigtAverage( values=[], weights=[] ):
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
	
def reussAverage( values=[], weights=[] ):
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
	
def vrhAverage( values=[], weights=[] ):
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
			K1	=	Bulk modulus of one phase of the mixture	(GPa)
			Mu1	=	Shear modulus of one phase of the mixture	(GPa)
			f1	=	Volume fraction of one phase of the mixture	(fraction, V/V)
			K2	=	Bulk modulus of the other phase of the mixture	(GPa)
			Mu2	=	Shear modulus of the other phase of the mixture	(GPa)
			f1	=	Volume fraction of one phase of the mixture	(fraction, V/V)
		
		Output Variable:
			Khs	=	Hashin Shtrikman bulk modulus			(GPa)
			
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
		
		A = 1 / (K2 - K1)
		B = f1 / (K1 + 4*Mu1/3)
		C = f2 / (A + B)
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
		
		A = 1 / (K2 - K1)
		B = f1 / (K1 + 4*Mu1/3)
		C = f2 / (A + B)
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
			K1	=	Bulk modulus of first phase of the mixture	(GPa)
			Mu1	=	Shear modulus of first phase of the mixture	(GPa)
			f1	=	Volume fraction of first phase of the mixture	(fraction, V/V)
			K2	=	Bulk modulus of the 2nd phase of the mixture	(GPa)
			Mu2	=	Shear modulus of the 2nd phase of the mixture	(GPa)
			f1	=	Volume fraction of second phase of the mixture	(fraction, V/V)
		
		Output Variable:
			Muhs	=	Hashin Shtrikman bulk modulus			(GPa)
			
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
		
		A = 1 / (Mu2 - Mu1)
		B1 = 2 * f1 * (K1 + 2*Mu1)
		B2 = 5 * Mu1 * (K1 + 4*Mu1/3)
		B = B1 / B2
		C = f2 / (A + B)

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
		
		A = 1 / (Mu2 - Mu1)
		B1 = 2 * f1 * (K1 + 2*Mu1)
		B2 = 5 * Mu1 * (K1 + 4*Mu1/3)
		B = B1 / B2
		C = f2 / (A + B)
		
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
	
def applyScalarFunctionToArrayVariables( function, parms=[] ):
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
		output = [function(parms[0][i],parms[1][i],
		parms[2][i]) for i in range(len(parms[0]))]
	elif numvars == 4:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i]) for i in range(len(parms[0]))]
	elif numvars == 5:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i],parms[4][i]) for i in range(len(parms[0]))]
	elif numvars == 6:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i],parms[4][i],
		arms[5][i]) for i in range(len(parms[0]))]
	elif numvars == 7:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i],parms[4][i],parms[5][i],
		parms[6][i]) for i in range(len(parms[0]))]
	elif numvars == 8:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i],parms[4][i],parms[5][i],
		parms[6][i],parms[7][i]) for i in range(len(parms[0]))]
	elif numvars == 9:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i],parms[4][i],parms[5][i],
		parms[6][i],parms[7][i],parms[8][i]) for i in range(len(parms[0]))]
	elif numvars == 10:
		output = [function(parms[0][i],parms[1][i],
		parms[2][i],parms[3][i],parms[4][i],parms[5][i],
		parms[6][i],parms[7][i],parms[8][i],parms[9][i]) for i in range(len(parms[0]))]
	else:
		mm.myStdErrorMessage( "Too many variables", 
		"Only up to 10 variables can be entered for a function" )
		mm.exitNow()
		
	return output
	
def matrix_K_calculator( constituentVolFrac = [], constituentK = [] ):
	"""
	The function computes the equivalent bulk modulus Kmin from the 
	contituent minerals within the matrix and their respective volume fractions
	using the Voigt-Reuss-Hill average.
	"""
	
	if len(constituentVolFrac) == len(constituentK):
	
		Kmin = vrhAverage(constituentK, constituentVolFrac)
		
	else:
		mm.myStdErrorMessage( "Length missmatch", 
		"The number of elements in weights do not match the numvber of constituents" )
		return
		
	return Kmin

def waterDensity( reservoirP, reservoirT ):
	"""
	This function computes the water density at reservoir conditions.
	The only input required are:
	
		reservoirP:	reservoir pressure		(MPa)
		reservoirT:	reservoir temperature		(degrees C)
		
	The output will be:
	
		waterDensity:	formation brine density		(g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP

	A1 = -80*T - 3.3*T*T + 0.00175*T*T*T + 489*P - 2*T*P
	A2 = 0.016*T*T*P - 1.3E-5*T*T*T*P - 0.333*P*P - 0.002*T*P*P
	Density = 1 + 1E-6*(A1 + A2)

	return Density
	
def brineDensity( salinity, reservoirP, reservoirT ):
	"""
	This function computes the brine density at reservoir conditions.
	The only input required are:
	
		salinity:	the salinity of formation brine	(weight fraction)
		reservoirP:	reservoir pressure		(MPa)
		reservoirT:	reservoir temperature		(degrees C)
		
	The output will be:
	
		brineDensity:	formation brine density		(g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	S = salinity
	
	densityWater = waterDensity( P, T )
	
	B1 = 80 + 3*T - 3300*S - 13*P + 47*P*S
	B2 = 300*P - 2400*P*S
	B3 = 0.668*S + 0.44*S*S

	brineDensity = densityWater + B3 + 1E-6*S*(B2 + T*B1)

	return brineDensity

def waterVelocity( reservoirP, reservoirT ):
	"""
	This function computes the compressional velocity of water at reservoir conditions.
	The only input required are:
	
		reservoirP:	reservoir pressure			(MPa)
		reservoirT:	reservoir temperature			(degrees C)
		
	The output will be:
	
		pVelWater:	formation brine velocity		(m/s)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP

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
	w[3][3] = -2.135E-6
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

	return pVelWater
	
	
def pVelBrine( salinity, reservoirP, reservoirT ):
	"""
	The function computes the compressional velocity of brine at reservoir conditions.
	The only input required are:
	
		salinity:	the salinity of formation brine	        (weight fraction)
		reservoirP:	reservoir pressure			(MPa)
		reservoirT:	reservoir temperature			(degrees C)
		
	The output will be:
	
		pVelBrine:	formation brine velocity		(m/s)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	S = salinity

	pVelWater = waterVelocity( P, T )
				
	A = 1170 - 9.6*T + 0.055*T*T - 8.5E-5*T*T*T + 2.6*P - 0.0029*T*P - 0.0476*P*P
	B = 780 - 10*P + 0.16*P*P
	pVelBrine = pVelWater + S*A + (S**1.5)*B - 1820*S*S
	
	return pVelBrine
	
	
def brineBulkModulus( salinity, reservoirP, reservoirT ):
	"""
	This function computes the bulk mudulus of brine at reservoir conditions.
	The only input required are:
	
		salinity:		the salinity of formation brine			(weight fraction)
		reservoirP:		reservoir pressure				(MPa)
		reservoirT:		reservoir temperature				(degrees C)
		
	The output will be:
	
		Kbrine:			Brine bulk modulus at reservoir conditions	(GPa)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	S = salinity

	
	Rbrine = brineDensity(S,P,T)
	Vbrine = pVelBrine(S,P,T)
	
	Kbrine = Rbrine * Vbrine*Vbrine * 1E-6
	
	return Kbrine
	
def gasDensity( gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the density of gas at reservoir conditions.
	The only input required are:
	
		gasSpecificGravity:	the specific gravity of gas			(API at 15.6C and 1 atmosphere)
		reservoirP:		reservoir pressure				(MPa)
		reservoirT:		reservoir temperature				(degrees C)
		
	The output will be:
	
	        Rgas:			gas density at reservoir conditions	        (g/cm^3)

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
	
		gasSpecificGravity:		the specific gravity of gas		(API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure			(MPa)
		reservoirT:			reservoir temperature			(degrees C)
		
	The output will be:
	
		Rgas:				gas density at reservoir conditions	(g/cm^3)

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

		oilApiGravity:			API gravity of oil			(API at 15.6C and 1 atmosphere)
		gasOilRatio:			gas to oil ratio GOR			(litre/litre)
		gasSpecificGravity:		ratio of gas density to air density	(API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure			(MPa)
		reservoirT:			reservoir temperature			(degrees C)
		
	The output will be:
	
		Roil:				oil density at reservoir conditions	(g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	Rg = gasOilRatio
	R0 = 141.5/( oilApiGravity + 131.5 )     # Convert from oil API gravity to specific gravity

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

		oilApiGravity:			API gravity of oil				(API at 15.6C and 1 atmosphere)
		gasOilRatio:			gas to oil ratio				(litre/litre)
		gasSpecificGravity:		ratio of gas density to air density	        (API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure				(MPa)
		reservoirT:			reservoir temperature				(degrees C)
		
	The output will be:
	
		Vp_oil:				oil p-wave velocity at reservoir conditions	(g/cm^3)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	Rg = gasOilRatio
	R0 = 141.5/( oilApiGravity + 131.5 )		# Convert from oil API gravity to specific gravity

	B1 = (2.4 * Rg * sqrt(G/R0) + T + 17.8)**1.175
#	B1 = (2.495 * Rg * sqrt(G/R0) + T + 17.8)**1.175
	B0 = 0.972 + 0.00038 * B1			# formation volume factor

	Rps = R0 / ((1 + 0.001 * Rg) * B0)
	Rs = (R0 + 0.0012 * Rg * G) / B0

	A = 2096*sqrt(Rps / (2.6 - Rps)) - 3.7*T + 4.64*P
	B = 0.0115*(sqrt((18.33/Rps) - 16.97) - 1)*T*P
	Vp_oil = A + B

	return Vp_oil

def oilBulkModulus( oilApiGravity, gasOilRatio, gasSpecificGravity, reservoirP, reservoirT ):
	"""
	This function computes the oil denisty at reservoir conditions.
	The only input required are:

		oilApiGravity:			API gravity of oil				(API at 15.6C and 1 atmosphere)
		gasOilRatio:			gas to oil ratio				(litre/litre)
		gasSpecificGravity:		ratio of gas density to air density		(API at 15.6C and 1 atmosphere)
		reservoirP:			reservoir pressure				(MPa)
		reservoirT:			reservoir temperature				(degrees C)
		
	The output will be:
	
		Koil:				oil bulk modulus at reservoir conditions	(Gpa)

	Kumar Geohorizons Jan 2006
	"""
	T = reservoirT
	P = reservoirP
	G = gasSpecificGravity
	Rg = gasOilRatio
	R0 = 141.5/( oilApiGravity + 131.5 )    # Convert from api gravity to specific gravity gm/cc

	Voil = oilPwaveVelocity(R0,Rg,G,P,T)    # Oil p-wave velocity m/s
	Roil = oilDensity(R0,Rg,G,P,T)          # Oll density g/cm^3

	Koil = Roil * Voil*Voil * 1E-6

	return Koil
	
def gassmannModeling( initWaterSat, finalWaterSat, 
		      initOilSat, finalOilSat, 
		      oilApiGravity, 
		      gasSpecificGravity, 
		      gasOilRatio, 
		      initReservoirP, finalReservoirP,
		      initEffectiveP,
		      initReservoirT, finalReservoirT, 
		      initWaterSalinity, finalWaterSalinity, 
		      initPorosity, finalPorosity,
		      volumeShale,
		      initVp, 
		      initVs, 
		      initDensity,
		      dryFrameKScale = 1.0,
		      dryFrameMuScale = 1.0,
		      alterDryFrameMod = 1,
		      nBiot = 1.0,
		      printVar = 0,
		      sigma = 1.0,
		      **kwargs):
	"""
	This function attempts to group all the variables in one place but still be general enough
	for testing of various reservoir conditions. This assumes there is an initial state described
	by the variables with the prefix 'init' and a final state defined by variables prefixed with 
	'final'. Currently some variables are fixed for the exercise and will not have the designated
	prefixes. It is assumed the fluid in the system is composed of gas, oil and brine at some
	saturation, such that, (water saturation)+(oil saturation)+(gas saturation) = 1. This means only two
	of these saturations need to be specified in both the inital and final states. I have chosen 
	use the water and oil saturations so gas saturation will be calculated with the previous mentioned
	relationship. The output of this exercise will be the final Vp,Vs and density for the specified
	final state. If you want the dry-frame moduli to remain unchanged through a change in pressure then
	set 'alterDryFramMod' to anything but 1. Setting 'alterDryFramMod' to 1 will mean you want
	to ajust the dry-frame moduli using Hertz-Mindlin theory.

	Setting printVar to 1 if you want all intermediate results printed, 2 for partial and
	anything else for no results printed. The default is no printing.
	
	The units and brief definition of input variables follows:
	
	
		(Variable Name)				(Brief Description)		       	(Units Expected)
		initWaterSat, finalWaterSat:		water saturation		       	(V/V volume fraction)
		initOilSat, finalOilSat:		oil saturation				(V/V volume fraction)
		oilApiGravity:				API gravity of oil			(API at 15.6C and 1 atmosphere)
		gasSpecificGravity:			Gas density / Air density		(API at 15.6C and 1 atmosphere)
		gasOilRatio:				gas to oil ratio			(litre/litre)
		initReservoirP, finalReservoirP:	Reservoir pore pressure			(MPa)
		initEffectiveP, finalEffectiveP:	Reservoir effective pressure		(MPa)
		initReservoirT, finalReservoirT:	Reservoir temperature			(degrees C)
		initWaterSalinity: 			Initial brine salinity			(ppm will require conversion)
		finalWaterSalinity:			Final brine salinity			(ppm will require conversion)
		initPorosity, finalPorosity:		Porosity				(V/V volume fraction)
		volumeShale				Fraction shale in matrix		(V/V volume fraction)
		initVp:					initial P-wave velocity			(m/s)
		initVs:					initial S-wave velocity			(m/s)
		initDensity:				initial bulk density			(g/cm^3)
		
	Some of these will require conversion to units excepted by other functions called, such as, oil API
	gravity will be converted to specific gravity (SG = 141.5/(API + 131.5)) and salinity will go from ppm
	to weight fraction (wf = ppm * 10^-6). The conversion will be done in this function.
	
	The ouput variables are:
	
		(Variable Name)				(Brief Descripton)			(Units)
		finalVp					final P-wave velocity			(m/s)
		finalVs					final S-wave velcoity			(m/s)
		finalDenisty				final bulk density			(g/cm^3)
	
	Currently This function assumes a sand-shale environment so certain fixed parameters are set below
	"""
	# I make the assumtion that Plith (litho-static pressure) is constant so Peff and Ppore must be related
	# through the expression Peff = Plith - n*Ppore, where: Peff is the effective pressure, Ppore is the 
	# reservoir pore pressure and n (Biot/Willis constant) which is often taken to be 1.0.
	# This can be set by the user. The default is n=1.0. This means we can compute the effective
	# pressure as:
	finalEffectiveP = initEffectiveP - nBiot*(finalReservoirP - initReservoirP)

	if printVar == 1 or printVar == 2:
		print "\n\nInput variables:\n" 
		print "     Initial water saturation:      "+str(initWaterSat)+"  (V/V volume fraction)"
		print "     Final water saturation:        "+str(finalWaterSat)+"  (V/V volume fraction)"
		print "     Initial oil saturation:        "+str(initOilSat)+"  (V/V volume fraction)"
		print "     Final oil saturation:          "+str(finalOilSat)+"  (V/V volume fraction)"
		print "     Oil API gravity:               "+str(oilApiGravity)+"  (API at 15.6C and 1 atmosphere)" 
		print "     Gas specific gravity:          "+str(gasSpecificGravity)+"  (API at 15.6C and 1 atmosphere)"
		print "     Gas to Oil ratio:              "+str(gasOilRatio)+"  (litre/litre)"
		print "     Initial pore pressure:         "+str(initReservoirP)+"  (MPa)"
		print "     Final pore pressure:           "+str(finalReservoirP)+"  (MPa)"
		print "     Initial effective pressure:    "+str(initEffectiveP)+"  (MPa)"
		print "     Final effective pressure:      "+str(finalEffectiveP)+"  (MPa)"
		print "     Initial reservoir temperature: "+str(initReservoirT)+"  (degress C)"
		print "     Final reservoir temperature:   "+str(finalReservoirT)+"  (degress C)"
		print "     Initial water salinity:        "+str(initWaterSalinity)+"  (ppm)"
		print "     Final water salinity:          "+str(finalWaterSalinity)+"  (ppm)"
		print "     Initial porosity:              "+str(initPorosity)+"  (V/V volume fraction)"
		print "     Final porosity:                "+str(finalPorosity)+"  (V/V volume fraction)"
		print "     Shale volume fraction:         "+str(volumeShale)+"  (V/V volume fraction)"
		print "     Initial P-wave velocity:       "+str(initVp)+"  (m/s) = "+str(initVp*3.280833)+" (ft/s)"
		print "     Initial S-wave velocity:       "+str(initVs)+"  (m/s) = "+str(initVs*3.280833)+" (ft/s)"
		print "     Initial bulk density:          "+str(initDensity)+"  (g/cc)"
		print "     Biot-Willis coefficient n:     "+str(nBiot)+"  (none)"
		print "     "
	
	# Fixed parameters for sand-shale environment (Mavko et al. 1998)
	
	Kclay = 20.9			# Clay bulk modulus (Gpa)
	Kquartz = 36.6			# Quartz bulk modulus (Gpa)
	MuClay = 6.67			# Clay shear modulus (Gpa)
	MuQuartz = 45.0			# Quartz shear modulus (Gpa)
	Rho_clay = 2.58			# Clay bulk density (g/cm^3)
	Rho_quartz = 2.65		# Quartz bulk density (g/cm^3)
	

	if printVar == 1 or printVar == 2:
		print "Fixed parameters:\n"
		print "     Clay bulk modulus:     "+str(Kclay)+"  (Gpa)"
		print "     Quartz bulk modulus:   "+str(Kquartz)+"  (Gpa)"
		print "     Clay shear modulus:    "+str(MuClay)+"  (Gpa)"
		print "     Quartz shear modulus:  "+str(MuQuartz)+"  (Gpa)"
		print "     Clay density:          "+str(Rho_clay)+"  (g/cc)"
		print "     Quartz density:        "+str(Rho_quartz)+"  (g/cc)"
		print "     "
	
	# Assumptions, conversions and computed quantities
	
	volumeClay = volumeShale * 0.7				# Assumption Vclay = 0.7 * Vshale in the matrix
	volumeQuartz = 1 - volumeClay				# Assume quartz makes up the rest of the matrix
											

	if initWaterSat + initOilSat > 1.0 or initWaterSat < 0.0 or initOilSat < 0.0 :
		print "Impossible initial fluid saturations input."
		return 0
	else:
		initGasSat = 1-initWaterSat-initOilSat		# Compute initial gas saturation

	if finalWaterSat + finalOilSat > 1.0 or finalWaterSat < 0.0 or finalOilSat < 0.0 :
		print "Impossible final fluid saturations input."
		return 0
	else:
		finalGasSat = 1-finalWaterSat-finalOilSat	# Compute final gas saturation
		
	initWaterSalinity = initWaterSalinity * 1E-6		# Convert salinity from PPM to weight fraction
	finalWaterSalinity = finalWaterSalinity * 1E-6		# Convert salinity from PPM to weight fraction
	oilSpecificGravity = 141.5/(oilApiGravity+131.5)	# Convert oil API gravity to oil specific gravity
	init_Mu_sat = get_Mu_from_Vs_Rho( initVs, initDensity )
	init_Mu_dry = init_Mu_sat                               # Gassmann assumption
	init_K_sat = get_K_from_Vp_Vs_Rho( initVp, initVs, initDensity )
	
	if printVar == 1 or printVar == 2:
		print "Assumptions and Conversions:\n"
		print "     Volume fraction clay:   "+str(volumeClay)+"  Assumed to be .7*Vsh (V/V)"
		print "     Volume fraction quartz: "+str(volumeQuartz)+"  Assumed to be 1-Vclay (V/V)"
		print "     Initial gas saturation: "+str(initGasSat)+"  (V/V)"
		print "     Final gas saturation:   "+str(finalGasSat)+"  (V/V)"
		print "     Initial water salinity: "+str(initWaterSalinity)+"  Converted from ppm to weight fraction (g/cc)"
		print "     Final water salinity:   "+str(finalWaterSalinity)+"  Converted from ppm to weight fraction (g/cc)"
		print "     Oil specific gravity:   "+str(oilSpecificGravity)+"  Convert oil API grav to oil specific gravity"
		print "     "
		print "Initial bulk elastic properties from (Vp,Vs,Rho):\n"
		print "     Initial shear modulus:  "+str(init_Mu_sat)+"  (Gpa)"
		print "     Initial bulk modulus:   "+str(init_K_sat)+"  (Gpa)"
		print "     "

	# Compute bulk matrix properties (note sand-clay-shale assumptions above)
	
	Kmatrix = matrix_K_calculator([volumeClay, volumeQuartz],[Kclay, Kquartz])
	MuMatrix = matrix_K_calculator([volumeClay, volumeQuartz],[MuClay, MuQuartz])
	RhoMatrix = voigtAverage([Rho_clay, Rho_quartz],[volumeClay, volumeQuartz])
	
	if printVar == 1:
		print "Matrix properties:\n"
		print "     Matrix bulk modulus :    "+str(Kmatrix)+"  (Gpa)"
		print "     Matrix shear modulus:    "+str(MuMatrix)+"  (Gpa)"
		print "     Matrix bulk density :    "+str(RhoMatrix)+"  (g/cc)"
		print "     "
		
	# Compute water/brine properties

	initWaterDensity = waterDensity( initReservoirP, initReservoirT )
	initWaterVelocity = waterVelocity( initReservoirP, initReservoirT )

	finalWaterDensity = waterDensity( finalReservoirP, finalReservoirT )
	finalWaterVelocity = waterVelocity( finalReservoirP, finalReservoirT )
	
	initRhoBrine = brineDensity( initWaterSalinity, initReservoirP, initReservoirT )
	init_K_Brine = brineBulkModulus( initWaterSalinity, initReservoirP, initReservoirT )
	
	finalRhoBrine = brineDensity( finalWaterSalinity, finalReservoirP, finalReservoirT )
	final_K_Brine = brineBulkModulus( finalWaterSalinity, finalReservoirP, finalReservoirT )
	
	if printVar == 1:
		print "Water/Brine properties:\n"
		print "     Initial water density:  "+str(initWaterDensity)+"  (g/cc)"
		print "     Initial water velocity: "+str(initWaterVelocity)+"  (m/s)"
		print "     Final water density:    "+str(finalWaterDensity)+"  (g/cc)"
		print "     Final water velocity:   "+str(finalWaterVelocity)+"  (m/s)"
		print "     Initial brine density:  "+str(initRhoBrine)+"  (g/cc)"
		print "     Initial brine bulk mod: "+str(init_K_Brine)+"  (Gpa)"
		print "     Final brine density:    "+str(finalRhoBrine)+"  (g/cc)"
		print "     Final brine bulk mod:   "+str(final_K_Brine)+"  (Gpa)"
		print "     "
		
	# Compute oil properties
	
	initRhoOil = oilDensity( oilApiGravity, gasOilRatio, gasSpecificGravity, initReservoirP, initReservoirT )
	init_K_Oil = oilBulkModulus( oilApiGravity, gasOilRatio, gasSpecificGravity, initReservoirP, initReservoirT )
	
	finalRhoOil = oilDensity( oilApiGravity, gasOilRatio, gasSpecificGravity, finalReservoirP, finalReservoirT )
	final_K_Oil = oilBulkModulus( oilApiGravity, gasOilRatio, gasSpecificGravity, finalReservoirP, finalReservoirT )

	if printVar == 1:
		print "Oil properties:\n"
		print "     Initial oil density:    "+str(initRhoOil)+"  (g/cc)"
		print "     Initial oil bulk mod:   "+str(init_K_Oil)+"  (Gpa)"
		print "     Final oil density:      "+str(finalRhoOil)+"  (g/cc)"
		print "     Final oil bulk mod:     "+str(final_K_Oil)+"  (Gpa)"
		print "     "
	
	# Compute gas properties
	
	initRhoGas = gasDensity( gasSpecificGravity, initReservoirP, initReservoirT )
	init_K_Gas = gasBulkModulus( gasSpecificGravity, initReservoirP, initReservoirT )
	
	finalRhoGas = gasDensity( gasSpecificGravity, finalReservoirP, finalReservoirT )
	final_K_Gas = gasBulkModulus( gasSpecificGravity, finalReservoirP, finalReservoirT )

	if printVar == 1:
		print "Gas properties:\n"
		print "     Initial gas density:    "+str(initRhoGas)+"  (g/cc)"
		print "     Initial gas bulk mod:   "+str(init_K_Gas)+"  (Gpa)"
		print "     Final gas density:      "+str(finalRhoGas)+"  (g/cc)"
		print "     Final gas bulk mod:     "+str(final_K_Gas)+"  (Gpa)"
		print "     "
	
	# Compute bulk fluid properties 
	
	initRhoFluid = voigtAverage( [initRhoBrine,initRhoOil,initRhoGas], [initWaterSat,initOilSat,initGasSat] )
	init_K_Fluid = reussAverage( [init_K_Brine,init_K_Oil,init_K_Gas], [initWaterSat,initOilSat,initGasSat] )
	
	finalRhoFluid = voigtAverage( [finalRhoBrine,finalRhoOil,finalRhoGas], [finalWaterSat,finalOilSat,finalGasSat] )
	final_K_Fluid = reussAverage( [final_K_Brine,final_K_Oil,final_K_Gas], [finalWaterSat,finalOilSat,finalGasSat] )

	if printVar == 1:
		print "Bulk fluid properties:\n"
		print "     Initial fluid density:  "+str(initRhoFluid)+"  (g/cc)"
		print "     Initial fluid bulk mod: "+str(init_K_Fluid)+"  (Gpa)"
		print "     Final fluid density:    "+str(finalRhoFluid)+"  (g/cc)"
		print "     Final fluid bulk mod:   "+str(final_K_Fluid)+"  (Gpa)"
		print "     "
		
	# Dry frame properties
	
	init_K_dry = get_Kdry_FromGassmannEquation( init_K_sat, init_K_Fluid, Kmatrix, initPorosity )
	
	if alterDryFrameMod != 1:
		
		final_K_dry = init_K_dry		# This assumption needs to be tested
		final_Mu_dry = init_Mu_dry
		
	else:
		# Compute dry-frame moduli for pressure change based on Hertz-Mindlin theory
		final_K_dry,final_Mu_dry = dryFrameModuliChangeWithPressure( finalEffectiveP,
									     initEffectiveP,
									     init_K_dry,
									     init_Mu_dry,
									     initPorosity,
									     Kmatrix,
									     MuMatrix,
									     scaleK = sigma,
									     scaleMu = sigma )

	if printVar == 1:
		print "Dry frame properties:\n"
		print "     Initial dry frame bulk mod: "+str(init_K_dry)+"  (Gpa)"
		print "     Final dry frame bulk mod:    "+str(final_K_dry)+"  (Gpa)"
		print "     Initial dry frame shear mod: "+str(init_Mu_dry)+"  (Gpa)"
		print "     Final dry frame shear mod:   "+str(final_Mu_dry)+"  (Gpa)"
		print "     "
		
	# Computed initial bulk density from user variables compared to initial input density

	init_Rho_sat = voigtAverage( [initRhoFluid,RhoMatrix], [initPorosity,1.0-initPorosity] )
	densityRatio = initDensity/init_Rho_sat
	# if you wish to apply this correction factor (density ratio)  to computed values so that original density
	# is a fixed value at initial conditions set the following variable to 1 anything else will mean
	# use computed values.
	applyDensityCorrection = 1       
	if printVar == 1:
		print "Computed vs. Initial Input Density:\n"
		print "     Initial saturated input denisty: "+str(initDensity)+"  (g/cc)"
		print "     Computed saturated density:      "+str(init_Rho_sat)+"  (g/cc)"
		if applyDensityCorrection == 1:
			print "     Use input density for computation."
		else:
			print "     Computed density used for computation."
		print "     "
	# Final bulk properties
	
	final_K_sat = get_Ksat_FromGassmannEquation( final_K_dry, final_K_Fluid, Kmatrix, finalPorosity )

	final_Rho_sat = voigtAverage( [finalRhoFluid,RhoMatrix], [finalPorosity,1.0-finalPorosity] )
	if applyDensityCorrection == 1:
		final_Rho_sat = densityRatio * final_Rho_sat
		
	final_Mu_sat = final_Mu_dry		# This is a Gassmann assumption
	
	finalVs = get_Vs_from_Mu_Rho( final_Mu_sat, final_Rho_sat )
	finalVp = get_Vp_from_K_Mu_Rho( final_K_sat, final_Mu_sat, final_Rho_sat )

	if printVar == 1 or printVar == 2:
		print "Final bulk properties:\n"
		print "     Final bulk density:    "+str(final_Rho_sat)+"  (g/cc)"
		print "     Final bulk modulus:    "+str(final_K_sat)+"  (Gpa)"
		print "     Final shear modulus:   "+str(final_Mu_sat)+"  (Gpa)"
		print "     Final P-wave velocity: "+str(finalVp)+"  (m/s) = "+str(finalVp*3.280833)+"  (ft/s)"
		print "     Final S-wave velocity: "+str(finalVs)+"  (m/s) = "+str(finalVs*3.280833)+"  (ft/s)"
		print "     "
	
	return finalVp, finalVs, final_Rho_sat

def elonTemperature( TVDss ):
	"""
	Elon temperature in centigrade (C) as a function of depth subsea (TVDss)
	From PVT review January 2008

	Input variable:

	    TVDss            True vertical depth subsea             (m)

	Output variable:

	    C                Temperature in degrees centigrade      (degree C)
	"""

	C = 5.7666E-2*TVDss + 14.38944

	return C

def elonPorePressure( TVDss ):
	"""
	Elon pressure in Mega-Pascals (Pp) as a function of depth subsea in meters (TVDss)
	From PVT review January 2008 (pre-production wells)

	Input variable:

	    TVDss            True vertical depth subsea             (m)

	Output variable:

	    Pp               Pore pressure      (MPa)
	"""

	Pp = (TVDss - 1100)*7.8407E-3 + 11.3675

	return Pp

def passList2Function( function, args ):
	function(*args)
	return

def delKeyFromDictionary( inputDictionary, itemName ):
	
	outDictionary = dict( inputDictionary )
	del outDictionary[ itemName ]
	return outDictionary

def newLithoStaticPressure( densityLog, TVDss, waterDepth = 0.0, seaFloorDensity = -999.25, null = -999.25 ):
	"""
	************** Work in progress *******************
	
	This function integrates the density log (input as a list in gm/cc) ignoring all
	null values. If the sea floor denisty value is set to "null" (which is the default case)
	it will be set to the first ono-null value of the density log, otherwise, a linear trend is
	assumed between the sea floor density and the first density value of the log. The TVDss
	list is assumed to have the same dimension as the densityLog list. Be careful, this
	function assumes TVDss is positive downwards from sealevel which is opposite from Petrel.
	The TVDss log is assumed to contain no "null" values or only at points where the density log
	is also a null value. The default "null" value is -999.25 and can be altered.
	"""

	P = []
	g = 9.80665            # gravitational acceleration at earth surface in mks units (m/s^2)
	Dold = 0.0
	Zold = 0.0
	Pold = 0.0
	first = 0
	Zs = TVDss[0]
	Ze = TVDss[-1]
	Zw = float(waterDepth)
	Ds = densityLog[0]
	
	i = 0
	while True:
		if densityLog[i] != null:
			if seaFloorDensity == null:
				D0 = densityLog[i]         # setting seafloor density to
			Zs = TVDss[i]                      # first non-null density log value
			Ds = densityLog[i]
			break
		else:
			i+=1
			
	P0 = (Zs-Zw)*g*(D0+Ds)/2.0
	
	if len(densityLog) == len(TVDss) :
		for i in range(len(densityLog)):
			Znew = TVDss[i]
			if densityLog[i] == null:
				D = null
			else:
				D = densityLog[i] * 1.0E3           # Convert g/cm^3 to kg/m^3
				
			if D == null:
				P.append(null)
				
			elif Znew < waterDepth:
				P.append(0.0)
				Zold = Znew
				
			elif Znew >= waterDepth and Znew <= Zs:
				Pnew = D0*g*(Znew-waterDepth) + (Ds-D0)*g*(Znew-waterDepth)**2/(2*(Zs-Zw))
				P.append(Pnew * 1.0E-6)             # Store pressure as MPa
				Pold = Pnew
				
			elif Znew > Zs and first == 0:
				first = 1
				Pold = (Zs-Zw)*g*(D0+Ds)/2.0
				Pnew = Pold + D*g*(Znew-Zold)       # Compute pressure from water bottom in Pa
				P.append(Pnew * 1.0E-6)             # Store pressure as MPa
				Dold = D
				Zold = Znew
				Pold = Pnew
				
			elif Znew > Zs and first == 1:
				Pnew = Pold + D*g*(Znew-Zold)       # Compute pressure in Pa
				P.append(Pnew * 1.0E-6)             # Store pressure as MPa
				Dold = D
				Zold = Znew
				Pold = Pnew
				
			else:
				Pnew = Pold + Dold*g*(Znew-Ze)
				P.append(Pnew * 1.0E-6)             # Store pressure as MPa
				Dold = D
				Zold = Znew
				Pold = Pnew
				
	else:
		print " "
		print "**** Error in: lithoStaticPressure ****"
		print "     The length of density log must equal the length of the TVDss values."
		print " "

	return P

def lithoStaticPressure( densityLog, TVDss, waterDepth = 0.0, null = -999.25 ):
	"""
	This function integrates the density log (input as a list in gm/cc) ignoring all
	null values. The first density reading is assumed to be carried to the waterDepth TVDss = waterDepth
	The TVDss list is assumed to have the same dimension as the densityLog list.
	Be careful, this function assumes TVDss is positive downwards from sealevel
	which is opposite from Petrel. The default "null" value is -999.25 and can be altered.
	"""

	P = []
	g = 9.80665            # gravitational acceleration at earth surface in mks units (m/s^2)
	Dold = 0.0
	Zold = 0.0
	Pold = 0.0
	first = 0
	
	if len(densityLog) == len(TVDss) :
		for i in range(len(densityLog)):
			Znew = TVDss[i]
			if densityLog[i] == null:
				D = null
			else:
				D = densityLog[i] * 1.0E3           # Convert g/cm^3 to kg/m^3
				
			if D == null:
				P.append(null)
			elif Znew < waterDepth:
				P.append(0.0)
				Zold = Znew
			elif Znew > waterDepth and first == 0:
				first = 1
				Zold = waterDepth
				Pnew = Pold + D*g*(Znew-Zold)       # Compute pressure from water bottom in Pa
				P.append(Pnew * 1.0E-6)             # Store pressure as MPa
				Dold = D
				Zold = Znew
				Pold = Pnew
			else:
				Pnew = Pold + D*g*(Znew-Zold)       # Compute pressure in Pa
				P.append(Pnew * 1.0E-6)             # Store pressure as MPa
				Dold = D
				Zold = Znew
				Pold = Pnew
	else:
		print " "
		print "**** Error in: lithoStaticPressure ****"
		print "     The length of density log must equal the length of the TVDss values."
		print " "

	return P



def akiRichardsRpp( Vp1, Vs1, Rho1, Vp2, Vs2, Rho2 ):
	"""
	This program calculates the Aki-Richards approximate coefficients
	for PP reflection between two elastic half saces, as found in (Mavko 2009).
	I will use the same definition of the coeffients exactly as defined there
	though there appears to be obvious simplification to be made by
	rearrangemnent of the terms. The properites suffixed with 1 is assumed
	to be on top and the one suffixed by 2 on the bottom.

	The output of the program will be the zero-offset or normal incidence reflection
	coefficent R0, the second term will be the gradient term R1 and the final term
	will be the curvature term R2

	Input
	Vp1             Compressional Velocity in Layer 1          (m/s)
	Vs1             Shear Velocity in Layer 1                  (m/s)
	Rho1            Density in layer 1                         (g/cc)
	Vp2             Compressional Velocity in Layer 2          (m/s)
	Vs2             Shear Velocity in Larye 2                  (m/s)
	Rho2            Density in Layer 2                         (g/cc)

	Output
	R0              reflection strength intercept
	                Zero offset reflection coefficent
	R1              reflection strength gradient
	                controls mid range reflection strength
	R2              reflection strength curvature
	                control near critical reflection strength

        It is useful to recall that when modeled with spherical wavefronts the larger are
	quite different from what is predicted by plane wave theory which is what this is
	based on (Geophysics, Vol.73, No.2, March-April 2008, P.C7-C12,"Experimental
	verification of spherical-wave effect on the AVO response and implications for
	three-term inversion").

        where: Rpp = R0 + R1*(sin(theta)**2) + R2*(tan(theta)**2 - sin(theta)**2)
	       theta: is aproximately the p-wave incident angle
	"""
	Rho1 = Rho1 * 1.0e3      # convert density from g/cc to kg/m^3
	Rho2 = Rho2 * 1.0e3      # convert density from g/cc to kg/m^3

	dVp = Vp2 - Vp1
	dVs = Vs2 - Vs1
	dR  = Rho2  - Rho1

	AVp = (Vp2 + Vp1)/2.0
	AVs = (Vs2 + Vs1)/2.0
	AR  = (Rho2 + Rho1)/2.0

	R0 = (dVp/AVp + dR/AR)/2.0
	R1 = dVp/(2.0*AVp) - (2.0*AVs**2/AVp**2)*(2.0*dVs/AVs + dR/AR)
	R2 = dVp/(2.0*AVp)

	return R0, R1, R2

def extractValueAtReference(originalValues, originalReference, desiredReference ):
	"""
	This function extracts a value at a user specified desired reference point
	from an original list of values and reference that may not contain a point
	at the desired reference (normally this would be something like a depth
	reference). The program will linearly interpolate a result to the desired point.
	All null values within the original values and reference will be removed before
	interpolation will be done. The input reference values are assumed to be
	monotonically increasing.

	"""
	originalValues = removeNull(originalReference, originalValues)
	originalReference = removeNull(originalValues,removeNull(originalReference))
	originalValues = removeNull(originalValues)

	if len(originalValues) != len(originalReference):
		print " "		
		print "**** Error from extractValueAtReference ****"
		print "     input reference arrays are not of the same length"
		print " "
	else:
		originalValues = removeNull(originalReference, originalValues)
		originalReference = removeNull(originalValues,removeNull(originalReference))
		originalValues = removeNull(originalValues)

		oldVal = 0.0
		oldRef = 0.0
		newVal = 0.0
		newRef = 0.0
		endIndex = len(originalReference) - 1
		
		if desiredReference < originalReference[0]:
				oldRef = originalReference[0]
				oldVal = originalValues[0]
				newRef = originalReference[1]
				newVal = originalValues[1]

				alpha = (desiredReference - oldRef)/(newRef - oldRef)
				desiredValue = alpha * (newVal - oldVal) + oldVal

				return desiredValue

		elif desiredReference > originalReference[endIndex]:
				oldRef = originalReference[endIndex-1]
				oldVal = originalValues[endIndex-1]
				newRef = originalReference[endIndex]
				newVal = originalValues[endIndex]

				alpha = (desiredReference - oldRef)/(newRef - oldRef)
				desiredValue = alpha * (newVal - oldVal) + oldVal

				return desiredValue
			
		for i in range(len(originalReference)):
			
			if originalReference[i] >= desiredReference:
				oldRef = originalReference[i-1]
				oldVal = originalValues[i-1]
				newRef = originalReference[i]
				newVal = originalValues[i]

				alpha = (desiredReference - oldRef)/(newRef - oldRef)
				desiredValue = alpha * (newVal - oldVal) + oldVal

				return desiredValue

		print " "
		print "**** Error in extractValueAtReference ****"
		print "     No value extracted and end of list reached."
		print "     Check to make sure your input original reference is monotonically increasing."
		print " "
		return
		
	
	
	
def removeNull( referenceList, inputList=['noInputSupplied'], null=-999.25):
	"""
	This funtion attempts to remove values in the input list where the referencList
	has a null value. If only one list is given then nulls of the referenceList
	will be removed.
	"""
	if inputList[0]=="noInputSupplied":
		inputList = referenceList
	
	outList = []
	for i in range(len(referenceList)):
		if referenceList[i] != null:
			outList.append(inputList[i])
			
	return outList

def medianFilt( input, filtLengthSamples=5, null=-999.25):
	"""
	This function runs a 'filtLengthSamples' point median file on the input
	list/array (currently only lists are considered). Null values will be
	returned as null values. The output will be the median filtered list/array.
	Hopefully you will have entered meaningful numbers in your input.

	"""
	listLength = len(input)
	A = []
	indexA = []
	B = []
	C = []
	
	output = [input[i] for i in range(listLength)]
	for i in range(listLength):
		if input[i] != null:
			A.append(input[i])
			indexA.append(i)

	lengthNoNulls = len(A)
	
	if lengthNoNulls < filtLengthSamples:
		print " "	
		print "***** Error Median Filter *****"
		print "      Input list/array does not have enough non-null elements to filter."
		print " "
		return
	
	if filtLengthSamples%2 == 0:
		midFilt = int(filtLengthSamples/2)
		oddEven = 0
	else:
		midFilt = int(filtLengthSamples/2) + 1
		oddEven = 1

	for i in range(filtLengthSamples):
		B.append(A[i])
		C.append(A[i])
		
	C.sort()
	
	if filtLengthSamples%2 == 0:
		med = float(C[midFilt-1]+C[midFilt])/2.0
		C = []
	else:
		med = float(C[midFilt-1])
		C = []

	for i in range(midFilt):
		A[i] = med
		
	del B[0]

	for i in range(midFilt,lengthNoNulls-midFilt):
		B.append(A[midFilt+i-oddEven])
		C = [B[j] for j in range(len(B))]
		C.sort()
		
		if filtLengthSamples%2 == 0:
			med = float(C[midFilt-1]+C[midFilt])/2.0
			C = []
		else:
			med = float(C[midFilt-1])
			C = []
		A[i] = med
		del B[0]

	for i in range(lengthNoNulls-midFilt,lengthNoNulls):
		A[i] = med

	for i in range(lengthNoNulls):
		output[indexA[i]] = A[i]

	return output

def despike( input, medFiltLength=5, deviationFromMedian=2.0, null=-999.25 ):
	"""
	This function removes spikes as determined by deviation from the median
	reference curve. The statistics are computed from the difference between
	the input and its median filtered version. Null values are ignored

	"""
	medInput = medianFilt( input, medFiltLength, null )
	sum = 0.0
	count = 0
	
	A = [input[i] for i in range(len(input))]
	
	diff = [A[i] - medInput[i] for i in range(len(input))]
	
	for i in range(len(A)):
		if input != null:
			sum += diff[i]*diff[i]
			count += 1

	stdDev = sqrt(sum)/float(count)

	startDelIndex = 0
	endDelIndex = 0
	startRefIndex = 0
	endRefIndex = 0
	delFlag = 0

	for i in range(len(A)):
		if A[i] != null:
			if fabs(diff[i]) >= fabs(medInput[i]) + deviationFromMedian*stdDev and delFlag == 0:
				startDelIndex = i
				startRefIndex = i - 1
				delFlag = 1
			elif fabs(diff[i]) < fabs(medInput[i]) + deviationFromMedian*stdDev and delFlag == 1:
				endDelIndex = i - 1
				endRefIndex = i
				delFlag = 0
				for j in range(startDelIndex,endDelIndex + 1):
					alpha = (endRefIndex - j) / (endRefIndex - startRefIndex)
					if A[j] != null:
						A[j] = A[startRefIndex]*alpha + (1.0-alpha)*A[endRefIndex]
					

	return A

	

	
def readOneRowFromAscii( fileName, rowNumber, delimiter=" " ):
    """
    Read a single row from an ascii file and output an list with
    each element of the array being delimiter separated strings from
    your specified row (default delimiter is spaces).

    """
    f = open(fileName,'r')
    count = 0
    for line in f:
	    count += 1
	    if count == rowNumber:
		    if delimiter == " ":
			    output = line.strip().split()
			    f.close()
			    return output
		    else:
			    output = line.strip().split(delimiter)
			    f.close()
			    return output

    print ""
    print "**** Error from readOneRowFromAscii ****"
    print ""
    print "      Row number was not reached your file maybe too short"
    print ""
    return	
		

			
def readRowsFromAscii( fileName, startRowNumber=1, endRowNumber=-999 , delimiter=" "):
    """
    Read specified rows from an ascii file and output an list with
    each element of the list being one of the defined columns. It is
    assuemed that all the rows have the same number of elements.
    Currently the defaults are to start at the first row and read
    to the end of file and the elements are spaces delimited. All
    these can be set by the user.

    """
    f = open(fileName,'r')
    
    count = 0
    first = 0
    
    for line in f:
	    
	    count += 1
	    
	    if endRowNumber == -999 and count >= startRowNumber :
		    
		    if first == 0:
			    output = []
			    if delimiter == " ":
				    tempList = line.strip().split()
				    for i in range(len(tempList)):
					    output.append([tempList[i]])
					    first = 1
			    else:
				    tempList = line.strip().split(delimiter)
				    for i in range(len(temp)):
					    output.append([tempList[i]])
					    first = 1
		    else:
			    if delimiter == " ":
				    tempList = line.strip().split()
				    for i in range(len(tempList)):
					    output[i].append(tempList[i])
			    else:
				    tempList = line.strip().split(delimiter)
				    for i in range(len(tempList)):
					    output[i].append(tempList[i])
					    
	    elif count >= startRowNumber and count <= endRowNumber:
		    
		    if first == 0:
			    output = []
			    if delimiter == " ":
				    tempList = line.strip().split()
				    for i in range(len(tempList)):
					    output.append([tempList[i]])
					    first = 1
			    else:
				    tempList = line.strip().split(delimiter)
				    for i in range(len(tempList)):
					    output.append([tempList[i]])
					    first = 1
		    else:
			    if delimiter == " ":
				    tempList = line.strip().split()
				    for i in range(len(tempList)):
					    output[i].append(tempList[i])
			    else:
				    tempList = line.strip().split(delimiter)
				    for i in range(len(tempList)):
					    output[i].append(tempList[i])

	    elif endRowNumber != -999 and  count > endRowNumber :
		    f.close()
		    return output
	    
    if endRowNumber != -999:
	    print " "
	    print "**** Warning from readRowsFromAscii ****"
	    print " "
	    print "      End of file reached (hopefully this is what you wanted."
	    print " "

	    f.close()
	    return output
    else:
	    f.close()
	    return output

def pLithE2LinearFit( TVDss ):
	"""
	A rough liner fit to the deeper values (800m to ~1100m TVDss) of Plith (MPa)
	and TVDss (m) to be used only for a rough pass for litho static pressure
	estimate.
	"""

	m = 0.022143           # slope
	b = -3.14021           # intercept
	
	pLith = m*TVDss + b    # linear approximation

	return pLith

			
def linearFit( x, y ):
	"""
	This function uses numpy routines perform a linear fit of two one
	dimensional numerical python lists or numpy arrays (x,y) of equal
	length. The output will be a list with the intercept and gradient
	of the fitted line.

	"""

	if type(x) == list and type(y) == list:
		
		lx = len(x)
		ly = len(y)
	
		if lx != ly:
			print " "
			print "**** Error in linearFitList ***"
			print "     input lists are not of the same length"
			print " "
			return
		else:
			X = np.asarray(x)
			Y = np.asarray(y)
			A = np.vstack([X,np.ones(lx)]).T
			W = np.linalg.lstsq(A,Y)[0]
			
	elif type(x) == np.ndarray and type(y) == np.ndarray:
		
		lx = len(x)
		ly = len(y)
	
		if lx != ly:
			print " "
			print "**** Error in linearFitList ***"
			print "     input numpy arrays are not of the same length"
			print " "
			return'Rho_Del_Sgas_SwSoil100'
		else:
			A = np.vstack([x,np.ones(lx)]).T
			W = np.linalg.lstsq(A,y)[0]
			
	w = list(W)
	return w

def findKey( inputDictionary, pattern ):
	count = 0
	for k in inputDictionary:
		if pattern in k:
			print "'"+str(k)+"'"
			count += 1
	if count == 0:
		print "**** No matches found for "+str(pattern)+" ****"
	else:
		print str(count)+" matches found"
	return

def findNotKey( inputDictionary, pattern ):
	count = 0
	for k in inputDictionary:
		
		if pattern in k:
			a=0
		else:
			print "'"+str(k)+"'"
			count += 1
	if count == 0:
		print "**** No miss-matches found for "+str(pattern)+" ****"
	else:
		print str(count)+" miss-matches found"
	return

			     
def fndPrnKeyVal( inputDictionary, pattern ):
	count = 0
	for k in inputDictionary:
		if pattern in k:
			val = inputDictionary[k]
			print "'"+str(k)+"' = "+str(val)
			count += 1
	if count == 0:
		print " "
		print "**** No matches found for "+str(pattern)+" ****"
	else:
		print " "
		print str(count)+" matches found"
	return

def fndPrnNotKeyVal( inputDictionary, pattern ):
	count = 0
	for k in inputDictionary:
		if pattern in k:
			a=0
		else:
			val = inputDictionary[k]
			print "'"+str(k)+"' = "+str(val)
			count += 1
	if count == 0:
		print " "
		print "**** No miss-matches found for "+str(pattern)+" ****"
	else:
		print " "
		print str(count)+" miss-matches found"
	return

def gman( finalWaterSat,
	  finalOilSat,
	  finalReservoirP,
	  finalReservoirT,
	  finalWaterSalinity,
	  finalPorosity,
	  initWaterSat = 0.1,
	  initOilSat = 0.9,
	  initReservoirP = 10.846,
	  initEffectiveP = 8.810,
	  initReservoirT = 73.99,
	  initWaterSalinity = 160000.0,
	  initPorosity = 0.28,
	  oilApiGravity = 31.3,
	  gasSpecificGravity = 0.939,
	  gasOilRatio =38.6, 
	  volumeShale = 0.06,
	  initVp = 2542, 
	  initVs = 1340, 
	  initDensity = 2.168,
	  dryFrameKScale = 1.0,
	  dryFrameMuScale = 1.0,
	  **kwargs):
	"""
	This function loops through a few allowed user chosen final parameter in Gassman substitution
	that is a list with other parameters as scalars. Currently if finalWaterSat is a list
	you must enter a list for finalOilSat and the sum of two correspoding elements of the
	list must be <= 1.0. This is because the remainder is assumed to be gas so that
	waterSat + oilSat + gasSat = 1.0

	
	"""
	if type(finalReservoirT) == list and type(finalWaterSat) == list and type(finalOilSat) == list and type(finalWaterSalinity) == list:
		
		if len(finalReservoirT) == len(finalWaterSat) and len(finalWaterSat) == len(finalOilSat) and len(finalOilSat) == len(finalWaterSalinity):
			Vp = []
			Vs = []
			R = []
			for i in range(len(finalReservoirT)):
				a,b,c = gassmannModeling( initWaterSat, finalWaterSat[i], 
							  initOilSat, finalOilSat[i], 
							  oilApiGravity, 
							  gasSpecificGravity, 
							  gasOilRatio, 
							  initReservoirP, finalReservoirP,
							  initEffectiveP,
							  initReservoirT, finalReservoirT[i], 
							  initWaterSalinity, finalWaterSalinity[i], 
							  initPorosity, finalPorosity,
							  volumeShale,
							  initVp, 
							  initVs, 
							  initDensity,
							  dryFrameKScale = 1.0,
							  dryFrameMuScale = 1.0 )
				Vp.append(a)
				Vs.append(b)
				R.append(c)
		else:
			print " "
			print "**** Error in gman ****"
			print "     If you want to vary temperature, salinity and saturation together,"
			print "     you must have the arrays be of equal length."
			print " "

			return
			


	
	elif type(finalWaterSat) == list and type(finalOilSat) == list:
		Vp = []
		Vs = []
		R = []
		for p in range(len(finalWaterSat)):
			a,b,c = gassmannModeling( initWaterSat, finalWaterSat[p], 
						  initOilSat, finalOilSat[p], 
						  oilApiGravity, 
						  gasSpecificGravity, 
						  gasOilRatio, 
						  initReservoirP, finalReservoirP,
						  initEffectiveP,
						  initReservoirT, finalReservoirT, 
						  initWaterSalinity, finalWaterSalinity, 
						  initPorosity, finalPorosity,
						  volumeShale,
						  initVp, 
						  initVs, 
						  initDensity,
						  dryFrameKScale = 1.0,
						  dryFrameMuScale = 1.0 )
			Vp.append(a)
			Vs.append(b)
			R.append(c)

	elif type(finalWaterSat) == list and type(finalOilSat) != list or (
		type(finalWaterSat) != list and type(finalOilSat) == list):

		print " "
		print "**** Error in gman ****"
		print "     Either your oil or water saturation is not a list."
		print "     Both must be either lists or scalars."
		print " "

	elif type(finalReservoirP) == list:
		Vp = []
		Vs = []
		R = []
		for p in finalReservoirP:
			a,b,c = gassmannModeling( initWaterSat, finalWaterSat, 
						  initOilSat, finalOilSat, 
						  oilApiGravity, 
						  gasSpecificGravity, 
						  gasOilRatio, 
						  initReservoirP, p,
						  initEffectiveP,
						  initReservoirT, finalReservoirT, 
						  initWaterSalinity, finalWaterSalinity, 
						  initPorosity, finalPorosity,
						  volumeShale,
						  initVp, 
						  initVs, 
						  initDensity,
						  dryFrameKScale = 1.0,
						  dryFrameMuScale = 1.0 )
			Vp.append(a)
			Vs.append(b)
			R.append(c)

	

	elif type(finalReservoirT) == list and type(finalWaterSat) != list and type(finalOilSat) != list:
		Vp = []
		Vs = []
		R = []
		for p in finalReservoirT:
			a,b,c = gassmannModeling( initWaterSat, finalWaterSat, 
						  initOilSat, finalOilSat, 
						  oilApiGravity,
						  gasSpecificGravity, 
						  gasOilRatio, 
						  initReservoirP, finalReservoirP,
						  initEffectiveP,
						  initReservoirT, p, 
						  initWaterSalinity, finalWaterSalinity, 
						  initPorosity, finalPorosity,
						  volumeShale,
						  initVp, 
						  initVs, 
						  initDensity,
						  dryFrameKScale = 1.0,
						  dryFrameMuScale = 1.0 )
			Vp.append(a)
			Vs.append(b)
			R.append(c)

	
				
	elif type(finalWaterSalinity) == list:
		Vp = []
		Vs = []
		R = []
		for p in finalWaterSalinity:
			a,b,c = gassmannModeling( initWaterSat, finalWaterSat, 
						  initOilSat, finalOilSat, 
						  oilApiGravity, 
						  gasSpecificGravity, 
						  gasOilRatio, 
						  initReservoirP, finalReservoirP,
						  initEffectiveP,
						  initReservoirT, finalReservoirT, 
						  initWaterSalinity, p, 
						  initPorosity, finalPorosity,
						  volumeShale,
						  initVp, 
						  initVs, 
						  initDensity,
						  dryFrameKScale = 1.0,
						  dryFrameMuScale = 1.0 )
			Vp.append(a)
			Vs.append(b)
			R.append(c)

	elif type(finalPorosity) == list:
		Vp = []
		Vs = []
		R = []
		for p in finalPorosity:
			a,b,c = gassmannModeling( initWaterSat, finalWaterSat, 
						  initOilSat, finalOilSat, 
						  oilApiGravity, 
						  gasSpecificGravity, 
						  gasOilRatio, 
						  initReservoirP, finalReservoirP,
						  initEffectiveP,
						  initReservoirT, finalReservoirT, 
						  initWaterSalinity, finalWaterSalinity, 
						  initPorosity, p,
						  volumeShale,
						  initVp, 
						  initVs, 
						  initDensity,
						  dryFrameKScale = 1.0,
						  dryFrameMuScale = 1.0 )
			Vp.append(a)
			Vs.append(b)
			R.append(c)

	else:
		Vp,Vs,R = gassmannModeling( initWaterSat, finalWaterSat, 
					  initOilSat, finalOilSat, 
					  oilApiGravity, 
					  gasSpecificGravity, 
					  gasOilRatio, 
					  initReservoirP, finalReservoirP,
					  initEffectiveP,
					  initReservoirT, finalReservoirT, 
					  initWaterSalinity, finalWaterSalinity, 
					  initPorosity, finalPorosity,
					  volumeShale,
					  initVp, 
					  initVs, 
					  initDensity,
					  dryFrameKScale = 1.0,
					  dryFrameMuScale = 1.0 )

	return Vp, Vs, R

def elon3rdOrderFitKdrytoPeff( inputPeff ):

	KdryCalc = -4.94741142e-04*inputPeff**3-5.54207296e-02*inputPeff**2+1.03781922e+00*inputPeff+2.14779425e-01

	return KdryCalc

def averageBands( bandProperty, numberBands=4, refProperty = [], null = -999.25):
	
	num = [0.0 for j in range(numberBands)]
	sum = [0.0 for j in range(numberBands)]
	band = []
	centre = []
	
	if refProperty == []:
		
		max = np.asarray(bandProperty).max()
		min = np.asarray(bandProperty).min()
		boundary = list(np.linspace(min,max,numberBands+1))
		for i in range(numberBands):
			centre.append(float((boundary[i]+boundary[i+1])/2.0))
			for j in range(len(bandProperty)):
				if bandProperty[j] >= boundary[i] and  bandProperty[j] < boundary[i+1]:
					sum[i] += bandProperty[j]
					num[i] += 1.0
	else:
		if len(bandProperty) == len(refProperty):
			
			max = np.asarray(refProperty).max()
			min = np.asarray(refProperty).min()
			print max,"   ",min
			boundary = list(np.linspace(min,max,numberBands+1))
			print boundary
			for i in range(numberBands):
				centre.append(float((boundary[i]+boundary[i+1])/2.0))
				for j in range(len(bandProperty)):
					if refProperty[j] >= boundary[i] and refProperty[j] < boundary[i+1]:
						sum[i] += bandProperty[j]
						num[i] += 1.0

		else:
			print " "
			print "*** Error in avarageBands ****"
			print "    Input reference list and property to be banded"
			print "    must be the same length."
			print " "

			return

	for i in range(numberBands):
		if num[i] == 0.0:
			band.append(null)
		else:
			band.append(sum[i]/num[i])

	return [band, centre]

def temp(a, inDict = {}, b = 2, c = 'hello',x=4,y=5,z=6, **kwargs ):
	if len(inDict) == 0:
		print " "
		print "No input dictionary."
		print " "
		print "a = "+str(a)+" , b = "+str(b)+" , c = "+str(c)		
	else:
		print 'inDict = ',inDict
	print "a = "+str(a)+" , b = "+str(b)+" , c = "+str(c)
	print "x = "+str(x)+" , z = "+str(y)+" , z = "+str(z)
	print " "
	print "Input alternate arguments follow:"
	for k in kwargs:
		print k+" = "+str(kwargs[k])
	print " "
	print "Arbitrary 3 other assigned variables a,b and c:"
       	print "a = "+str(a)+" , b = "+str(b)+" , c = "+str(c)
	print locals()
	return locals()

def HM_scale(Ppore,s=1.0):
	gmanVars = {'alterDryFrameMod': 1,
	'dryFrameScale': 1.0,
	'finalOilSat': 0.9,
	'finalPorosity': 0.28,
	'finalReservoirP': 10.846,
	'finalReservoirT': 73.988551311764709,
	'finalWaterSalinity': 160000.0,
	'finalWaterSat': 0.1,
	'gasOilRatio': 38.6,
	'gasSpecificGravity': 0.939,
	'initDensity': 2.168,
	'initEffectiveP': 8.81,
	'initOilSat': 0.9,
	'initPorosity': 0.28,
	'initReservoirP': 10.846,
	'initReservoirT': 73.988551311764709,
	'initVp': 2542.0,
	'initVs': 1340.0,
	'initWaterSalinity': 160000.0,
	'initWaterSat': 0.1,
	'nBiot': 1.0,
	'dryFrameKScale': 1.0,
	'dryFrameMuScale': 1.0,
	'oilApiGravity': 31.3,
	'printVar': 0,
	'volumeShale': 0.06}
	gmanVars['dryFrameKScale'] = s
	gmanVars['dryFrameMuScale'] = s
	gmanVars['finalReservoirP'] = Ppore
	
	Vp,Vs,Rho = gassmannModeling(**gmanVars)
	return Vp

def Mu_ModHertzMind_scale(Peff, Scale):
	Pnew = Peff
	Pref = 8.81
	Kref = 4.613721530809661
	MuRef = 4.012
	PhiRef = 0.28
	Kmatrix = 35.7
	MuMatrix = 39.8
	scaleK = Scale
	scaleMu = Scale

	if isinstance(Peff,(int, long, float, complex)):
		Pnew = Peff
		a,b = dryFrameModuliChangeWithPressure( Pnew, Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, scaleK )
		K = a
		Mu = b
		
	elif type(Peff) == np.ndarray:
		K = np.array([])
		Mu = np.array([])
		for i in range(len(Peff)):
			a,b = dryFrameModuliChangeWithPressure( Pnew[i], Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, scaleK )
			K=np.append(K,a)
			Mu=np.append(Mu,b)
			
	elif type(Peff) == list:
		K = []
		Mu = []
		for i in range(len(Peff)):
			a,b = dryFrameModuliChangeWithPressure( Pnew[i], Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, scaleK )
			K.append(a)
			Mu.append(b)
	else:
		mm.myStdErrorMessage( "Error: K_Mu_HM_scale ", (
			"The input variables are not consistent for this function" ))
		return
			

	return Mu

def K_ModHertzMind_scale(Peff, Scale):
	Pnew = Peff
	Pref = 8.81
	Kref = 4.613721530809661
	MuRef = 4.012
	PhiRef = 0.28
	Kmatrix = 35.7
	MuMatrix = 39.8
	scaleK = Scale
	scaleMu = Scale

	if isinstance(Peff,(int, long, float, complex)):
		Pnew = Peff
		a,b = dryFrameModuliChangeWithPressure( Pnew, Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, scaleK )
		K = a
		Mu = b
		
	elif type(Peff) == np.ndarray:
		K = np.array([])
		Mu = np.array([])
		for i in range(len(Peff)):
			a,b = dryFrameModuliChangeWithPressure( Pnew[i], Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, scaleK )
			K=np.append(K,a)
			Mu=np.append(Mu,b)
			
	elif type(Peff) == list:
		K = []
		Mu = []
		for i in range(len(Peff)):
			a,b = dryFrameModuliChangeWithPressure( Pnew[i], Pref, Kref, MuRef, PhiRef, Kmatrix, MuMatrix, scaleK )
			K.append(a)
			Mu.append(b)
	else:
		mm.myStdErrorMessage( "Error: K_Mu_HM_scale ", (
			"The input variables are not consistent for this function" ))
		return
			

	return K

def Vp_ModHertzMind_scale(Ppore, Scale):

	gmanVars = {'finalWaterSat' : 0.1,
		    'finalOilSat' : 0.9,
		    'finalReservoirP' : 8.810,
		    'finalReservoirT' : 73.99,
		    'finalWaterSalinity' : 160000.0,
		    'finalPorosity' : 0.28,
		    'initWaterSat' : 0.1,
		    'initOilSat' : 0.9,
		    'initReservoirP' : 10.846,
		    'initEffectiveP' : 8.810,
		    'initReservoirT' : 73.99,
		    'initWaterSalinity' : 160000.0,
		    'initPorosity' : 0.28,
		    'oilApiGravity' : 31.3,
		    'gasSpecificGravity' : 0.939,
		    'gasOilRatio' : 38.6,
		    'volumeShale' : 0.06,
		    'initVp' : 2542,
		    'initVs' : 1340,
		    'initDensity' : 2.168,
		    'dryFrameKScale' : 1.0,
		    'dryFrameMuScale' : 1.0,
		    'dryFrameScalar' : 1.0,
		    'sigma' : 1.0}

	gmanVars['sigma'] = float(Scale)
	
	if isinstance(Ppore,(int, long, float, complex)):
		gmanVars['finalReservoirP'] = float(Ppore)
		a,b,c = gassmannModeling( **gmanVars )
		Vp = a
		Vs = b
		Rho = c
		
	elif type(Ppore) == np.ndarray:
		Vp = np.zeros(len(Ppore))
		Vs = np.zeros(len(Ppore))
		Rho = np.zeros(len(Ppore))
		i = 0
		for p in Ppore:
			gmanVars['finalReservoirP'] = float(p)
			a,b,c = gassmannModeling( **gmanVars )
			Vp[i] = a
			Vs[i] = b
			Rho[i] = c
			i += 1
			
	elif type(Ppore) == list:
		Vp = []
		Vs = []
		Rho = []
		for p in Ppore:
			gmanVars['finalReservoirP'] = float(p)
			a,b,c = gassmannModeling( **gmanVars )
			Vp.append(a)
			Vs.append(b)
			Rho.append(c)
	else:
		print " "
		print "**** Error in Vp_ModHertzMind_scale: ****"
		print "     The input variables are not consistent for this function"
		print " "
		return
			

	return Vp
	
def Vs_ModHertzMind_scale(Ppore, Scale):
	gmanVars = {'finalWaterSat' : 0.1,
		    'finalOilSat' : 0.9,
		    'finalReservoirP' : 8.810,
		    'finalReservoirT' : 73.99,
		    'finalWaterSalinity' : 160000.0,
		    'finalPorosity' : 0.28,
		    'initWaterSat' : 0.1,
		    'initOilSat' : 0.9,
		    'initReservoirP' : 10.846,
		    'initEffectiveP' : 8.810,
		    'initReservoirT' : 73.99,
		    'initWaterSalinity' : 160000.0,
		    'initPorosity' : 0.28,
		    'oilApiGravity' : 31.3,
		    'gasSpecificGravity' : 0.939,
		    'gasOilRatio' : 38.6,
		    'volumeShale' : 0.06,
		    'initVp' : 2542,
		    'initVs' : 1340,
		    'initDensity' : 2.168,
		    'dryFrameKScale' : 1.0,
		    'dryFrameMuScale' : 1.0,
		    'dryFrameScalar' : 1.0,
		    'sigma' : 1.0}

	gmanVars['sigma'] = float(Scale)
	
	if isinstance(Ppore,(int, long, float, complex)):
		gmanVars['finalReservoirP'] = float(Ppore)
		a,b,c = gassmannModeling( **gmanVars )
		Vp = a
		Vs = b
		Rho = c
		
	elif type(Ppore) == np.ndarray:
		Vp = np.zeros(len(Ppore))
		Vs = np.zeros(len(Ppore))
		Rho = np.zeros(len(Ppore))
		i = 0
		for p in Ppore:
			gmanVars['finalReservoirP'] = float(p)
			a,b,c = gassmannModeling( **gmanVars )
			Vp[i] = a
			Vs[i] = b
			Rho[i] = c
			i += 1
			
	elif type(Ppore) == list:
		Vp = []
		Vs = []
		Rho = []
		for i in Ppore:
			gmanVars['finalReservoirP'] = float(p)
			a,b,c = gassmannModeling( **gmanVars )
			Vp.append(a)
			Vs.append(b)
			Rho.append(c)
	else:
		print " "
		print "**** Error in Vp_ModHertzMind_scale: ****"
		print "     The input variables are not consistent for this function"
		print " "
		return
			

	return Vs

def Rho_ModHertzMind_scale(Ppore, Scale):
	gmanVars = {'finalWaterSat' : 0.1,
		    'finalOilSat' : 0.9,
		    'finalReservoirP' : 8.810,
		    'finalReservoirT' : 73.99,
		    'finalWaterSalinity' : 160000.0,
		    'finalPorosity' : 0.28,
		    'initWaterSat' : 0.1,
		    'initOilSat' : 0.9,
		    'initReservoirP' : 10.846,
		    'initEffectiveP' : 8.810,
		    'initReservoirT' : 73.99,
		    'initWaterSalinity' : 160000.0,
		    'initPorosity' : 0.28,
		    'oilApiGravity' : 31.3,
		    'gasSpecificGravity' : 0.939,
		    'gasOilRatio' : 38.6,
		    'volumeShale' : 0.06,
		    'initVp' : 2542,
		    'initVs' : 1340,
		    'initDensity' : 2.168,
		    'dryFrameKScale' : 1.0,
		    'dryFrameMuScale' : 1.0,
		    'dryFrameScalar' : 1.0,
		    'sigma' : 1.0}

	gmanVars['sigma'] = float(Scale)
	
	if isinstance(Ppore,(int, long, float, complex)):
		gmanVars['finalReservoirP'] = float(Ppore)
		a,b,c = gassmannModeling( **gmanVars )
		Vp = a
		Vs = b
		Rho = c
		
	elif type(Ppore) == np.ndarray:
		Vp = np.zeros(len(Ppore))
		Vs = np.zeros(len(Ppore))
		Rho = np.zeros(len(Ppore))
		i = 0
		for p in Ppore:
			gmanVars['finalReservoirP'] = float(p)
			a,b,c = gassmannModeling( **gmanVars )
			Vp[i] = a
			Vs[i] = b
			Rho[i] = c
			i += 1
			
	elif type(Ppore) == list:
		Vp = []
		Vs = []
		Rho = []
		for i in Ppore:
			gmanVars['finalReservoirP'] = float(p)
			a,b,c = gassmannModeling( **gmanVars )
			Vp.append(a)
			Vs.append(b)
			Rho.append(c)
	else:
		print " "
		print "**** Error in Vp_ModHertzMind_scale: ****"
		print "     The input variables are not consistent for this function"
		print " "
		return
			

	return Rho

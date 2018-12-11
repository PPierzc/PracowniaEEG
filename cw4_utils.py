import numpy as np
import mtmvar

def coeffs(signal, low, high):
	return [ mtmvar.mult_AR(signal, order, 1) for order in range(low, high + 1) ]

def akaike(p, V, shape):
	return np.log(np.linalg.det(V)) + 2 * p * shape[0] ** 2 / shape[1]

def get_aic_func(signal, low, high):
	cfs = coeffs(signal, low, high)
	aic = [akaike(p + 1, cf, signal.shape) for p, cf in enumerate(cfs)]
	return np.array(aic)

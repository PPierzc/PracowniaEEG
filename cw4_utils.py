import numpy as np
import PracowniaEEG.mtmvar as mtmvar

def coeffs(signal, low, high):
  return [ mtmvar.mult_AR(signal, order, 1) for order in range(low, high + 1) ]

def get_aic_func(signal, low, high):
  cfs = coeffs(signal, low, high)
  aic = [akaike(len(params), covars, signal.shape) for params, covars in cfs]
  return aic

def akaike(p, V, shape):
  """
  p: order
  V: covars matrix
  shape: input shape
  """

  return np.log(np.linalg.det(V)) + 2 * p * shape[0] * shape[0] / shape[1]

def best_akaike(signal, low, high):
  akaike_scores = get_aic_func(signal, 1, 20)
  best_index = akaike_scores.index(min(akaike_scores))
  order = low + best_index
  params, covars = mtmvar.mult_AR(signal, order, 1)
  return order, params, covars


def z_transform(z, f):
	T = 0.1
	F = sum([f[k] * z**(-k) for k in range(len(f)])
	return F


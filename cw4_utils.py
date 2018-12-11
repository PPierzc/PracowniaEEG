import numpy as np

\mathrm{AIC}(p)=\mathrm{ln}(\det({V}))+2\frac{pk^2}{N}

def akaike(p, V, shape):
	return np.log(np.linalg.det(V)) + 2 * p * shape[0] ** 2 / shape[1]


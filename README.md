# Complex-Curve-Fitting
Method from E.C. Levy's Paper: Complex curve fitting algorithm by E.C.Levy in IRE Transactions on Automatic Control AC-4 (1959)
## Example Usage:
Below snippet shows how to use the function on example data. The angle ```phase``` has to be in degree, the frequency ```omega``` has to be in radiant per second.
```python
import numpy as np
from cplx_curve_fitting import *
#system
mag = np.array([1, 1, 1.02, 1.12, 1.24, 1.44, 2.27, 4.44, 8.17, 10.05, 5.56, 2.55, 1.45, 1.00])
phase = np.array([0, 5, 10, 24, 31, 39, 51.5, 50.5, 28, -6, -59, -76, -82, -84])
omega = np.array([0.0, 0.1, 0.2, 0.5, 0.7, 1.0, 2.0, 4.0, 7.0, 10, 20, 40, 70, 100])
# function call:
N,D = cplx_curve_fit(mag, phase, omega, 2, 2)
print("\n")
print(N)
print("\n")
print(D)
```

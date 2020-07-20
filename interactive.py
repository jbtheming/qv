import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytz
import math

url = 'https://raw.githubusercontent.com/jbtheming/qv/master/data/TYU20_20200617.csv'
df = pd.read_csv(url, index_col='quote_datetime', parse_dates=True, infer_datetime_format=True)
df.index = df.index.tz_localize(pytz.utc)
df = df.sort_index(ascending=True)

data_type='price'

# Tuning parameters
k = 30 # Size of the window over which to perform local averaging
w = .47
λ = 1/12
a0 = 12

if(data_type == 'price'):
    Z = df.apply(np.log)
else:
    Z = df

# Calculate the returns
changes = Z['bid'].diff().dropna()
s0 = changes.std(ddof=1) # ddof= Delta Degrees of Freedom. The divisor used in calculations is N - ddof
dt = 0.00069444440305233 # del is a reserved word in Python. Some quantty representing the time between observations. What scale should these units be?
u = k * dt # Length of time window, in units of dt
v0 = a0 * s0 * u ** w # Initial truncation parameter

changes_trunc = np.clip(changes, a_min=-v0, a_max=v0) # Why did we only clip the positive changes in original?
s = changes_trunc.std(ddof=1)
v = a0 * s * u ** w
theta = k * math.sqrt(dt)

# Locally averaged increments
weights = np.arange(1, k+1, dtype=float) / k
weights = np.minimum(weights, 1-weights) # Localizing - needs to satisfy f(0) = 0, f(1) = 0 and integral from 0, 1 of f^2 > 0
zBar = changes.rolling(k-1).apply(lambda x: np.dot(x, weights[:-1])).dropna()

# Locally averaged squared increments
squaredWeights = np.square(np.diff(weights))
zHat = changes.rolling(k-1).apply(lambda x: np.dot(x ** 2, squaredWeights)).dropna()

qvNoisy = changes.pow(2).sum()
qvStep = zBar.pow(2) - (0.5 * zHat)
qv = qvStep.sum() / (k*λ)
qvContinuous = qvStep.loc[zBar.abs() <= v].sum() / (k*λ)
iGamma = zHat.loc[zBar.abs() <= v].sum() * ((theta**2)/(2*k*λ))




# def calculate_quadvar(data, data_type='price'):
#     """ Calculate realized quadratic variation of a series of prices or rates.

# %calculate total quadratic variation QV, the continous part QVc and
# %integrated noise variance Igamma
#     """
    
# #df = df.reindex(pd.date_range(start=df.index.min(), end=df.index.max(), freq='1T'))      
# #df = df.interpolate(method='linear')     

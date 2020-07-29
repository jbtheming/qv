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
dt = 0.00069444440305233 # del is a reserved word in Python. Some quantity representing the time between observations. What scale should these units be?
u = k * dt # Length of time window, in units of dt
v0 = a0 * s0 * u ** w # Initial truncation parameter

changes_trunc = np.clip(changes, a_min=-v0, a_max=v0)
s = changes_trunc.std(ddof=1)
v = a0 * s * u ** w
theta = k * math.sqrt(dt)

# Locally averaged increments
weights = np.arange(1, k+1, dtype=float) / k
weights = np.minimum(weights, 1-weights) # Localizing - needs to satisfy f(0) = 0, f(1) = 0 and integral from 0, 1 of f^2 > 0
zBar = changes.rolling(k-1).apply(lambda window: np.dot(window, weights[:-1])).dropna()

# Locally averaged squared increments
squaredWeights = np.square(np.diff(weights)) # Not nan because Numpy, not Pandas
zHat = changes.rolling(k-1).apply(lambda window: np.dot(window ** 2, squaredWeights)).dropna()

qvNoisy = changes.pow(2).sum()
qvStep = zBar.pow(2) - (0.5 * zHat)
qv = qvStep.sum() / (k*λ)
# If the locally averaged price change is greater than 
# the standard deviation of the continuous part, then the price change is probably a jump
qvContinuous = qvStep.loc[zBar.abs() <= v].sum() / (k*λ)
iGamma = zHat.loc[zBar.abs() <= v].sum() * ((theta**2)/(2*k*λ))




# def calculate_quadvar(data, data_type='price'):
#     """ Calculate realized quadratic variation of a series of prices or rates.

# %calculate total quadratic variation QV, the continous part QVc and
# %integrated noise variance Igamma
#     """
    
# #df = df.reindex(pd.date_range(start=df.index.min(), end=df.index.max(), freq='1T'))      
# #df = df.interpolate(method='linear')     

"""
        /// <summary>
        /// Normalises an accumalated quadratic variation of a bond future's price to basis point per day volatility
        /// </summary>
        /// <param name="interval">The time interval over which the <paramref name="qv"/>qv</param> has accumalated.
        /// <param name="qv">The accumalated quadratic variation.</param>
        /// <param name="ctdDuration">The modified duration of the future's cheapest to deliver bond.</param>
        /// <returns></returns>
        public static double Normalise(TimeSpan interval, double qv, double ctdDuration) // Specific to bond futures, not quad var in general, so does not belong in here
        {
            if(qv < 0.0)
            {
                return 0.0;
            }
            // The proportion of a day that is in our accumalation interval
            double intervalProportion = interval.TotalSeconds / TimeSpan.FromDays(1).TotalSeconds;
            // The rate at which quad variation is accumulating, per unit time
            double accumalationRate = qv / intervalProportion;
            // Take square root and multiply by duration of the CTD bond for future in question
            double normalisedQv = Math.Sqrt(accumalationRate) * ctdDuration;
            // Mutltiply by 100 to get a basis point per day volatility
            return normalisedQv * 100.0;
        }


    public class TuningParameters
    {
        public static readonly TuningParameters Default = new TuningParameters();

        public double w { get; set; } = .47;
        public double Lambda { get; set; } = 1.0 / 12.0;
        public double a0 { get; set; } = 1.5;
        public int MaxObservations { get; set; } = 600;
        public double Theta { get; set; } = 30.0 / Math.Sqrt(60.0 * 24.0);
    }


        /// <summary>
        /// Return the difference between two dates, expressed as the decimal proportion of a day, and with single-second precision.
        /// e.g. 24 hours = 1.0, 1 hour = 1 / 24, 1 second = 1 / (60 * 60 * 24)
        /// </summary>
        /// <param name="startTime"></param>
        /// <param name="endTime"></param>
        /// <returns></returns>
        public static double CalculateTimeDelta(DateTimeOffset startTime, DateTimeOffset endTime)
        {
            return Math.Sqrt((endTime - startTime).TotalSeconds / TimeSpan.FromDays(1).TotalSeconds);
        }
"""
import numpy as np
import pandas as pd
import pytz

url = 'https://raw.githubusercontent.com/jbtheming/qv/master/data/TYU20_20200617.csv'
df = pd.read_csv(url, index_col='quote_datetime', parse_dates=True, infer_datetime_format=True)
df.index = df.index.tz_localize(pytz.utc)
df

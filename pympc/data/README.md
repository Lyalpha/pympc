# Data

Data files used by `pympc`

## obs_codes.npy

These are observatory codes used by the IAU, published by the Minor Planet Center 
[here](https://minorplanetcenter.net/iau/lists/ObsCodes.html), last retrieved on Mar 04 2022.

They were packaged as a numpy binary file for speed of access by `pympc`, using the following
recipe:

```python
# Assume the observatory codes data is already stored as an ascii file at /tmp/obscodes.txt
import numpy as np
from astropy.table import Table

tbl = Table.read("/tmp/obscodes.txt", format="ascii.fixed_width", col_starts=(0,4,13,21), col_ends=(3,12,20,29))
arr = np.array(tbl)
np.save("pympc/data/obs_codes.npy", arr, allow_pickle=False)
```
# adsorption
A python package for estimating the best isotherm (i.e., Langmuir, Freundlich, Temkin, Toth, Sips, DR) and kinetic models (i.e., PFO, PSO, Weber-Morris, Avrami, Bangham, Elovich) as well as Arrhenius parameters (i.e., Ea, A).

[Install](#Install) / [Usage](#Usage)

# Install
First, make sure you have a Python 3 environment installed.
To install from github:
```Python
pip install -e git+https://github.com/amvro23/adsorption/#egg=adsorption
```
Note: It might be useful to write "git+https://github.com/amvro23/adsorption/#egg=adsorption" if installing directly from a Python interpreter as # can be interpreted as a comment.

# Usage

```Python
from adsorption import (Isotherms, Kinetics, ModifiedArrhenius)
import numpy as np
import pandas as pd
```


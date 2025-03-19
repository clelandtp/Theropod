import pandas as pd
from bokeh.io import output_notebook, export_svg, export_png
from bokeh.plotting import figure, show, output_file, save
from bokeh.models import ColumnDataSource, HoverTool, NumeralTickFormatter, Label
from bokeh.palettes import Category10
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
import os
import glob
import numpy as np
import pandas as pd
import sys

path = os.getcwd()
outputname = sys.argv[1]
stori_files = glob.glob(os.path.join(path, "*_stori.pkl"))
print(stori_files)
os.mkdir('output')

slopes = []

for f in stori_files:
    df = pd.read_pickle(f, compression = 'xz')
    td = 300
    rv = 0.99

    df['SlopeComb2']=np.where((df['P Death1']-0>= td) & (df['SciLin R-squared1']>= rv), df['SciLin Slope1'], 
                                np.where((df['P Death2']-df['P Death1']>= td) & (df['SciLin R-squared2']>= rv), df['SciLin Slope2'],
                                         np.where((df['P Death3']-df['P Death2']>= td) & (df['SciLin R-squared3']>= rv), df['SciLin Slope3'],
                                                  np.where((df['P Death4']-df['P Death3']>= td) & (df['SciLin R-squared4']>= rv), df['SciLin Slope4'],
                                                           np.where((df['P Death5']-df['P Death4']>= td) & (df['SciLin R-squared5']>= rv), df['SciLin Slope5'],0)))))
    df = df[df.SlopeComb2>0]

    df['SlopeChoice']=np.where((df['P Death1']-0>= td) & (df['SciLin R-squared1']>= rv), "Slope1", 
                                np.where((df['P Death2']-df['P Death1']>= td) & (df['SciLin R-squared2']>= rv), "Slope2",
                                         np.where((df['P Death3']-df['P Death2']>= td) & (df['SciLin R-squared3']>= rv), "Slope3",
                                                  np.where((df['P Death4']-df['P Death3']>= td) & (df['SciLin R-squared4']>= rv), "Slope4",
                                                           np.where((df['P Death5']-df['P Death4']>= td) & (df['SciLin R-squared5']>= rv), "Slope5",0)))))

    df2 = df[['Frequency','m/z', 'Intensity', 'Scan', 'Ions', 'SlopeComb2', 'SlopeChoice']].copy()

    slopes.append(df2)


df = pd.concat(slopes)
df = df.sort_values('m/z')
df.to_csv("output/"+outputname+".csv",index=False)
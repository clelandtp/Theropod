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

def kde_sklearn(x, x_grid, bandwidth=0.001, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

def create_q(width=1500, height=500,
            main_title='Protein Mass Spectrum'):
    tooltips = [
        ('Mass','@Mass{0.0000}'),
        ('Density','@Density')
        ]
    q = figure(
        width=width, height=height,
        title = main_title,
        tools = 'xwheel_zoom,xpan,box_zoom,undo,reset,save',
        tooltips=tooltips,
        output_backend= "canvas",
        )
    return q
q = create_q()

def create_r(width=1500, height=500,
            main_title='Protein Mass Spectrum'):
    tooltips = [
        ('Mass','@Mass{0.0000}'),
        ('Density','@Density')
        ]
    r = figure(
        width=width, height=height,
        title = main_title,
        tools = 'xwheel_zoom,xpan,box_zoom,undo,reset,save',
        tooltips=tooltips,
        output_backend= "svg",
        )
    return r
r = create_r()

path = os.getcwd()
outputname = sys.argv[1]
stori_files = glob.glob(os.path.join(path, "*_stori.pkl"))
print(stori_files)
os.mkdir('output')

#storidfs = (pd.read_pickle(f, compression = 'xz') for f in stori_files)
#print(storidfs)
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
    #print(df)
    df2 = df[['Frequency','m/z', 'Intensity', 'Scan', 'Ions', 'SlopeComb2', 'SlopeChoice']].copy()

    slopes.append(df2)
#print(slopes)

df = pd.concat(slopes)
df = df.sort_values('m/z')
df.to_csv("output/"+outputname+".csv",index=False)

df = df[((df['m/z']<=637.29) | 
         (df['m/z']>=637.3))] 

df = df[((df['m/z']<=638.29) | 
         (df['m/z']>=638.297))] 

df = df[((df['m/z']<=650.29) | 
         (df['m/z']>=638.297))] 

df = df[((df['m/z']<=654.25) | 
         (df['m/z']>=655.25))] 

df = df[((df['m/z']<=1629) | 
         (df['m/z']>=1632.7))] 

#df = df[((df['m/z']<=1637) | 
 #        (df['m/z']>=1639.7))] 

#df = df[((df['m/z']<=1642.5) | 
 #       (df['m/z']>=1645))] 

df = df[((df['m/z']<=535.25) | 
         (df['m/z']>=535.35))] 

df = df[((df['m/z']<=663.4) | 
         (df['m/z']>=663.46))] 

df = df[((df['m/z']<=664.4) | 
         (df['m/z']>=664.46))] 

df = df[((df['m/z']<=616.1) | 
         (df['m/z']>=616.2))] 

df = df[((df['m/z']<=617.1) | 
        (df['m/z']>=617.2))]

#df = df[((df['m/z']<=1586.4) | 
 #        (df['m/z']>=1600))]

df = df[((df['m/z']<=637.2) | 
         (df['m/z']>=638.32))]
		 
ppm = 1.5
thresh = (ppm/1e6)*df['m/z']
df['mzcluster'] = df['m/z'].diff().gt(thresh).cumsum().add(1)

df = df.groupby("mzcluster").filter(lambda x: len(x) >= 2)

sns.set_style(rc = {'axes.facecolor': 'white'})

plt.figure()
sns.scatterplot(data = df, x = 'm/z', y = 'SlopeComb2')
plt.savefig('output/mzvslope.png')

plt.figure()
sns.histplot(data = df, x = 'SlopeChoice', binwidth = 1)
plt.savefig('output/slopechoice.png')

mzcluslabelsmean = df.groupby("mzcluster")['SlopeComb2'].mean().reset_index()
mzcluslabelsmean = mzcluslabelsmean.rename(columns={'SlopeComb2': 'SlopeComb2mzc'})

df = df.merge(mzcluslabelsmean, 'left')

df['Chargemzclust'] = np.round(df['SlopeComb2mzc'] /2.292161, decimals = 0)
df['Massmzclust'] = df['m/z']*df['Chargemzclust'] - df['Chargemzclust'] * 1.0078250319

mzclusterrange = df.groupby('mzcluster').agg(mzmin=('m/z', 'min'), mzmax=('m/z', 'max'))

ztest = 100
isotest = 10
zneigh = 10
zchange = np.arange(-ztest,ztest+1,1)
isofind = np.arange(-isotest,isotest+1,1)
zneigh2 = np.arange(-zneigh,zneigh+1,1)
cmzm = {(z,n):[] for z in zchange for n in zneigh2}
cmzm2 = {(z,iso):[] for z in zchange for iso in isofind}
cmzmneigh = {(z,n):[] for z in zchange for n in zneigh2}
neighpos = {(z,n):[] for z in zchange for n in zneigh2}
cmzmlabels = {(z,n):[] for z in zchange for n in zneigh2}
cmzmiso = {(z,iso):[] for z in zchange for iso in isofind} 
isopos = {(z,iso):[] for z in zchange for iso in isofind}
cmzmlabels2 = {(z,iso):[] for z in zchange for iso in isofind}
neighsum ={(z):[] for z in zchange}
isosum = {(z):[] for z in zchange}
totalsum = {(z):[] for z in zchange}

mzidx = pd.IntervalIndex.from_arrays(mzclusterrange['mzmin'],
                                 mzclusterrange['mzmax'], 
                                 closed="both")

for z in zchange:
    #for n in zneigh2:
     #   print("z", z,"neighbor", n)
      #  cmzm[z,n] = df['Chargemzclust']+z
       # cmzmneigh[z,n] = np.where(cmzm[z,n]>=5, (df['m/z']*cmzm[z,n]+(1.0078250319*n))/(cmzm[z,n]+n), 0)
        #neighpos[z,n] = mzidx.get_indexer(cmzmneigh[z,n])
        #cmzmlabels[z,n] =  np.where(neighpos[z,n] != -1, 1.3, 0)

    for iso in isofind:
        print("z", z, "isotopologue", iso)
        cmzm2[z,iso] = df['Chargemzclust']+z
        cmzmiso[z,iso] = np.where(cmzm2[z,iso]>=5, df['m/z'] + (1.0078250319*(iso/cmzm2[z,iso])), 0)
        isopos[z,iso] = mzidx.get_indexer(cmzmiso[z,iso])
        cmzmlabels2[z,iso] =  np.where(isopos[z,iso] != -1, 1, 0)


for z in zchange: 
    #neighsum[z] = sum([cmzmlabels[i] for i in cmzmlabels.keys() if i[0]==z])
    isosum[z] = sum([cmzmlabels2[i] for i in cmzmlabels2.keys() if i[0]==z])
    #totalsum= {k: neighsum.get(k, 0) + isosum.get(k, 0) for k in set(neighsum) | set(isosum)}
    
dftotsum = pd.DataFrame(isosum)
dftotsum = dftotsum.replace(0, 1.0)
dftotsum2 = pd.merge(df[['m/z', 'Chargemzclust']], dftotsum, left_index=True, right_index=True)
dftotsum2['cmzCheck'] = dftotsum2.loc[:, 0:].sum(axis = 1)
dftotsum2['cmzChargeadjust'] = dftotsum2.iloc[:, 2:-1].idxmax(axis=1).astype('int64')
dftotsum2['cmzChargeadjust2'] = np.where(dftotsum2['cmzCheck']==(2*np.max(ztest)+1), 0, dftotsum2['cmzChargeadjust'])
#dftotsum2['cmzfinalcharge'] = np.where(dftotsum2['cmzCheck']==(2*np.max(ztest)+1), dftotsum2['Chargemzclust']*0, dftotsum2['Chargemzclust']+dftotsum2['cmzChargeadjust2'])
dftotsum2['cmzfinalcharge'] = dftotsum2['Chargemzclust']+dftotsum2['cmzChargeadjust2']
dftotsum2['mzcluster'] = df['mzcluster']
df['cmzfinalcharge'] = dftotsum2['cmzfinalcharge']
df['Massmzclust2'] = df['m/z']*df['cmzfinalcharge'] - df['cmzfinalcharge'] * 1.0078250319
df.to_csv("output/"+outputname+"_final.csv",index=False)

mz = df['m/z'].to_numpy()
#mczkde_range = np.arange(x.min(), x.max(),0.1)
mzkde_range = np.arange(500, 2000,0.01)
mzpdf = kde_sklearn(mz, mzkde_range, bandwidth=0.01, rtol=1E-4)
df252 = pd.DataFrame({'m/z': mzkde_range, 'Intensity': mzpdf}, columns=['m/z', 'Intensity'])
df252['Intensity'] = df252['Intensity']
df252['relIntensity'] = (df252['Intensity']/(df252['Intensity'].max()))*100

sns.set_style(rc = {'axes.facecolor': 'white'})
plt.figure()
plt.plot(mzkde_range, mzpdf*10000)
plt.xlim(650,1200)
plt.savefig('output/mzspectrum.png')

plt.figure()
sns.histplot(dftotsum2, x = 'm/z', y = 'cmzfinalcharge', binwidth = [5,1])
plt.ylim(0,)
plt.savefig('output/mzvfinalcharge.png')

plt.figure()
sns.histplot(dftotsum2, x = 'cmzfinalcharge', binwidth = 1)
plt.xlim(0,)
plt.savefig('output/finalchargehistogram.png')

plt.figure()
sns.histplot(data=df, x= 'cmzfinalcharge', y = 'Massmzclust2',fill=True,  binwidth = [1,1000], cbar = True)
plt.savefig('output/finalchargevmass.png')


x = df['Massmzclust2'].to_numpy()
#mczkde_range = np.arange(x.min(), x.max(),0.1)
mczkde_range = np.arange(5000, 75000,0.05)
pdf = kde_sklearn(x, mczkde_range, bandwidth=0.05, rtol=1E-4)
df25 = pd.DataFrame({'Mass': mczkde_range, 'Density': pdf}, columns=['Mass', 'Density'])
df25['Density'] = df25['Density']*100000000
df25['relDensity'] = (df25['Density']/(df25['Density'].max()))*100

output_file(filename="output/plotcanvas1.html", title="Mass spectrum")
cds = ColumnDataSource(data=df25)

q = create_q()
maxIntens = df25['relDensity'].max()
#Main line
q.line(
    'Mass', 'relDensity',
    source = cds,
    color = 'black',# alpha = 0.8,
    line_width = 2
    )

#Format axis labels
def add_axis_labels(q):
    q.xaxis.axis_label = 'Mass'
    q.xaxis.axis_label_text_font_size = '10pt'
    q.xaxis.major_label_text_font_size = '9pt'

    q.yaxis.axis_label = 'Relative Density'
    q.yaxis.axis_label_text_font_size = '10pt'
    q.yaxis.major_label_text_font_size = '9pt'
    q.yaxis.formatter = NumeralTickFormatter(format='0.')
add_axis_labels(q)

show(q)
save(q)

df26 = df25[['Mass', 'Density']].copy()
df27 = df25[['Mass', 'Density']].copy()
df26['Mass'] = df26['Mass']+1.0078250319
df26 = df26.rename(columns={'Mass': 'mz', 'Density': 'intensity'})
df27 = df27.rename(columns={'Mass': 'mz', 'Density': 'intensity'})
df26.to_csv("output/"+outputname+'-MH.csv', index = False)
df27.to_csv("output/"+outputname+'-M.csv', index = False)
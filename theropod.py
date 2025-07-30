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
import argparse
import subprocess

def kde_sklearn(x, x_grid, bandwidth=0.001, **kwargs):
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)



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

parser = argparse.ArgumentParser()

parser.add_argument("filename", type=str, help="Add a filename")
parser.add_argument("--charge", type=int, default = 60, help="Set the number of charges to evaluate; Default = 60")
parser.add_argument("--neighbor", type=int, default = 10, help="Set the number of either charge or isotope neighbors; Default = 10")
parser.add_argument("--tod", type=int, default = 300, help="Set the minimum ion life time; Default = 300")
parser.add_argument("--rsquare", type=int, default = 0.99, help="Set the r-squared value; Default = 0.99")
parser.add_argument("--slopecal", type=int, default = 2.292161, help="Set the KDE slope calibration value; Default = 2.292161")
parser.add_argument("--minmass", type=int, default = 5000, help="Set the minimum mass for the final neutral mass spectrum; Default = 5000")
parser.add_argument("--maxmass", type=int, default = 75000, help="Set the maximum mass for the final neutral mass spectrum; Default = 75000")
parser.add_argument("--chargecorrecttype", type=int, default = 1, help="Set the charge correction type 1 = isotopes, 2 = charge neighbors; Default = 1")
parser.add_argument("--spectrumtype", type=int, default = 1, help="Set the spectrum type 1 = kernal density estimation, 2 = mass histogram; Default = 1")


args = parser.parse_args()

path = os.getcwd()
outputname = args.filename
stori_files = glob.glob(os.path.join(path, "*_stori.pkl"))
print(stori_files)
os.mkdir('output')

#storidfs = (pd.read_pickle(f, compression = 'xz') for f in stori_files)
#print(storidfs)
slopes = []

##Merge and filter pickle files

for f in stori_files:
    df = pd.read_pickle(f, compression = 'xz')
    td = args.tod
    rv = args.rsquare

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



## Perform charge calculation and refinement
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

df['Chargemzclust'] = np.round(df['SlopeComb2mzc'] /args.slopecal, decimals = 0)
df['Massmzclust'] = df['m/z']*df['Chargemzclust'] - df['Chargemzclust'] * 1.0078250319

mzclusterrange = df.groupby('mzcluster').agg(mzmin=('m/z', 'min'), mzmax=('m/z', 'max'))

ztest = args.charge
isotest = args.neighbor
zneigh = args.neighbor
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

if args.chargecorrecttype == 2:
    for z in zchange:
        for n in zneigh2:
            print("z", z,"neighbor", n)
            cmzm[z,n] = df['Chargemzclust']+z
            cmzmneigh[z,n] = np.where(cmzm[z,n]>=5, (df['m/z']*cmzm[z,n]+(1.0078250319*n))/(cmzm[z,n]+n), 0)
            neighpos[z,n] = mzidx.get_indexer(cmzmneigh[z,n])
            cmzmlabels[z,n] =  np.where(neighpos[z,n] != -1, 1, 0)

elif args.chargecorrecttype == 1:
    for z in zchange:
        for iso in isofind:
            print("z", z, "isotopologue", iso)
            cmzm2[z,iso] = df['Chargemzclust']+z
            cmzmiso[z,iso] = np.where(cmzm2[z,iso]>=5, df['m/z'] + (1.0078250319*(iso/cmzm2[z,iso])), 0)
            isopos[z,iso] = mzidx.get_indexer(cmzmiso[z,iso])
            cmzmlabels2[z,iso] =  np.where(isopos[z,iso] != -1, 1, 0)

if args.chargecorrecttype == 2:
    for z in zchange: 
        neighsum[z] = sum([cmzmlabels[i] for i in cmzmlabels.keys() if i[0]==z])
        
elif args.chargecorrecttype == 1:
    for z in zchange:
        isosum[z] = sum([cmzmlabels2[i] for i in cmzmlabels2.keys() if i[0]==z])
    #totalsum= {k: neighsum.get(k, 0) + isosum.get(k, 0) for k in set(neighsum) | set(isosum)}
    
if args.chargecorrecttype == 2:    
    dftotsum = pd.DataFrame(neighsum)

elif args.chargecorrecttype == 1:
    dftotsum = pd.DataFrame(isosum)

#dftotsum = dftotsum.replace(0, 1.0)
dftotsum2 = pd.merge(df[['m/z', 'Chargemzclust']], dftotsum, left_index=True, right_index=True)
dftotsum2['cmzCheck'] = dftotsum2.loc[:, 0:].sum(axis = 1)
dftotsum2['cmzChargeadjust'] = dftotsum2.iloc[:, 2:-1].idxmax(axis=1).astype('int64')
dftotsum2['cmzChargeadjust2'] = np.where(dftotsum2['cmzCheck']==(2*np.max(ztest)+1), 0, dftotsum2['cmzChargeadjust'])
#dftotsum2['cmzfinalcharge'] = np.where(dftotsum2['cmzCheck']==(2*np.max(ztest)+1), dftotsum2['Chargemzclust']*0, dftotsum2['Chargemzclust']+dftotsum2['cmzChargeadjust2'])
dftotsum2['cmzfinalcharge'] = dftotsum2['Chargemzclust']+dftotsum2['cmzChargeadjust2']
dftotsum2['mzcluster'] = df['mzcluster']
df['cmzfinalcharge'] = dftotsum2['cmzfinalcharge']
df['Massmzclust2'] = df['m/z']*df['cmzfinalcharge'] - df['cmzfinalcharge'] * 1.0078250319
df['ioncount']=abs((df['Chargemzclust']/df['cmzfinalcharge'])).replace([np.inf, -np.inf], 1).apply(np.ceil).astype(int)
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

plt.figure()
y=sns.histplot(df['ioncount'], binwidth = 1, discrete = True)
y.bar_label(y.containers[1])
plt.savefig('output/ioncount.png')



#mczkde_range = np.arange(x.min(), x.max(),0.1)
if args.spectrumtype == 1:
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
    df2 = df.loc[df.index.repeat(df.ioncount)]
    x = df2['Massmzclust2'].to_numpy()
    mczkde_range = np.arange(args.minmass, args.maxmass,0.05)
    pdf = kde_sklearn(x, mczkde_range, bandwidth=0.05, rtol=1E-4)
    df25 = pd.DataFrame({'Mass': mczkde_range, 'Density': pdf}, columns=['Mass', 'Density'])
    df25['Density'] = df25['Density']*100000000
    df25['relDensity'] = (df25['Density']/(df25['Density'].max()))*100
    output_file(filename="output/plotcanvasKDE.html", title="Mass spectrum")
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
        
elif args.spectrumtype == 2:
    def create_q(width=1500, height=500,
                main_title='Protein Mass Spectrum'):
        tooltips = [
            ('Mass','@Mass{0.0000}'),
            ('Count','@Count')
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
    x = df['Massmzclust2'].to_numpy()
    mcz_range = np.arange(args.minmass, args.maxmass,0.05)
    histo = np.histogram(x, mcz_range, weights = df['ioncount'])
    histo2 = histo[0]
    mcz_range2 = histo[1]
    mcz_range2 = mcz_range[0:-1]
    df25 = pd.DataFrame({'Mass': mcz_range2, 'Count': histo2}, columns=['Mass', 'Count'])
    df25['relCount'] = (df25['Count']/(df25['Count'].max()))*100
    output_file(filename="output/plotcanvasHistogram.html", title="Mass spectrum")
    cds = ColumnDataSource(data=df25)
    q = create_q()
    maxIntens = df25['relCount'].max()
    #Main line
    q.line(
        'Mass', 'relCount',
        source = cds,
        color = 'black',# alpha = 0.8,
        line_width = 2
        )

    #Format axis labels
    def add_axis_labels(q):
        q.xaxis.axis_label = 'Mass'
        q.xaxis.axis_label_text_font_size = '10pt'
        q.xaxis.major_label_text_font_size = '9pt'

        q.yaxis.axis_label = 'Relative Count'
        q.yaxis.axis_label_text_font_size = '10pt'
        q.yaxis.major_label_text_font_size = '9pt'
        q.yaxis.formatter = NumeralTickFormatter(format='0.')
    add_axis_labels(q)

    show(q)
    save(q)
    df26 = df25[['Mass', 'Count']].copy()
    df27 = df25[['Mass', 'Count']].copy()
    df26['Mass'] = df26['Mass']+1.0078250319
    df26 = df26.rename(columns={'Mass': 'mz', 'Count': 'intensity'})
    df27 = df27.rename(columns={'Mass': 'mz', 'Count': 'intensity'})
    df26.to_csv("output/"+outputname+'-MH.csv', index = False)
    df27.to_csv("output/"+outputname+'-M.csv', index = False)



script_path = os.path.join(os.path.dirname(__file__), "theropod_mzML.R")


subprocess.run([r"C:/PROGRA~1/R/R-4.3.1/bin/x64/Rscript.exe", script_path])
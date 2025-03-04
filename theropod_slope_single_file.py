
import time as tm
import os
import glob
import numpy as np
import pandas as pd
import scipy.signal
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.fft import rfft
from sklearn.tree import DecisionTreeRegressor
import warnings
import sys
warnings.filterwarnings('ignore')

start = tm.perf_counter()

conversionB = 148242866.439
calval = np.sqrt(conversionB)
calval

path = os.getcwd()
#pkl_files = glob.glob(os.path.join(path, "*.ftr"))

pkl_files = sys.argv[1]

print(pkl_files)

data = pd.read_feather(pkl_files)
#print(f)
data['Channel A'] = (data['Channel A']/32767)*5000
signal = data['Channel A'].to_numpy()
maxtime = (len(signal)-1)*500
time = (np.linspace(0, maxtime, len(signal)))/1000000
signal = signal
time = time
Ts = np.mean(np.diff(time)) #Sampling Interval
Fs = 1/Ts                   #Sampling Frequency
Fn = Fs/2                   #Nyquist Frequency
fs = 1/np.mean(np.diff(time))
hh = np.hanning(len(signal))
hh[0:384000]=[1]
signal22 = hh*signal
signal2 = np.pad(signal22, (0, 4*len(signal)), 'empty')
xdft = rfft(signal2)
xdft2 = np.abs(xdft)
frequency = np.linspace(0, Fn, len(xdft2))
mz = np.square(calval/frequency[1:,])
mzidx = np.where(np.logical_and(mz>=500,mz<=2000))
mz2 = mz[mzidx]
xdft3 = xdft2[mzidx]
frequency2 = frequency[mzidx]
mzdiff = np.mean(np.diff(mz2))
spectrumoutput = np.transpose(np.vstack((mz2, xdft3)))
cols=['m/z', 'intensity']
spectrum = pd.DataFrame(spectrumoutput, columns = cols)
spectrum = spectrum.sort_values('m/z')
xdft3hat2 = interp1d(mz2, xdft3, kind = 'quadratic')
xdft3hat3 = interp1d(frequency2, xdft3, kind = 'quadratic')
mznew = np.linspace(np.max(mz2),np.min(mz2),15000000)
xdftnew = xdft3hat2(mznew)
peak_idx, _ = find_peaks(np.abs(xdftnew), prominence = 15000, height = 15000, distance = None)
freqnew = calval/np.sqrt(mznew[peak_idx])
print(os.path.splitext(os.path.basename(pkl_files))[0], np.size(mznew[peak_idx]))
pslopes1 = []
pslopes2 = []
pslopes3 = []
pslopes4 = []
pslopes5 = []
pdeath1 = []
pdeath2 = []
pdeath3 = []
pdeath4 = []
pdeath5 = []
scislopeslin1 = []
scislopeslinr1 = []
scislopeslinp1 = []
scislopeslin2 = []
scislopeslinr2 = []
scislopeslinp2 = []
scislopeslin3 = []
scislopeslinr3 = []
scislopeslinp3 = []
scislopeslin4 = []
scislopeslinr4 = []
scislopeslinp4 = []
scislopeslin5 = []
scislopeslinr5 = []
scislopeslinp5 = []

for i in range (0, len(freqnew)):
    ss = freqnew[i:i+1]
    print(ss)
    realstori1 = signal * np.cos((-2*np.pi)*ss*time)
    realstori1 = np.transpose(realstori1)
    realstori1 = np.add.accumulate(realstori1)
    imagstori1 = -signal * np.sin((-2*np.pi)*ss*time)
    imagstori1 = np.transpose(imagstori1)
    imagstori1 = np.add.accumulate(imagstori1)
    realstori21 = np.square(realstori1)
    imagstori21 = np.square(imagstori1)
    magstori1 = np.sqrt(realstori21+imagstori21)
    magstori2 = magstori1.flatten()
    timehat1 = scipy.signal.decimate(time,10)
    timehat1 = scipy.signal.decimate(timehat1,10)
    maghat1 = scipy.signal.decimate(magstori2, 10)
    maghat1 = scipy.signal.decimate(maghat1, 10)
    xs = timehat1
    ys = maghat1
    n_seg = 5
    dys = np.gradient(ys, xs)
    rgr = DecisionTreeRegressor(max_leaf_nodes=n_seg, min_impurity_decrease = 400, min_samples_leaf = 2)
    rgr.fit(xs.reshape(-1, 1), dys.reshape(-1, 1))
    dys_dt = rgr.predict(xs.reshape(-1, 1)).flatten()
    _, indices=np.unique(dys_dt, return_index = True)
    b = dys_dt[np.sort(indices)]
    print(b)
    TOD = xs[np.sort(indices)]
    TOD = np.append(TOD, max(xs))
    print(TOD)
    indices2 = np.sort(indices)
    print(len(indices2))
    print(indices2)
    if len(indices2) == 1:
        pslopes1.append(b[0])
        pdeath1.append(TOD[1])
        maglinear1 = scipy.stats.linregress(xs[indices2[0]:len(xs)], ys[indices2[0]:len(ys)])
        scislopeslin1.append(maglinear1.slope)
        scislopeslinr1.append(maglinear1.rvalue**2)
        scislopeslinp1.append(maglinear1.pvalue)
        pslopes2.append(0)
        pdeath2.append(0)
        scislopeslin2.append(0)
        scislopeslinr2.append(0)
        scislopeslinp2.append(0)
        pslopes3.append(0)
        pdeath3.append(0)
        scislopeslin3.append(0)
        scislopeslinr3.append(0)
        scislopeslinp3.append(0)
        pslopes4.append(0)
        pdeath4.append(0)
        scislopeslin4.append(0)
        scislopeslinr4.append(0)
        scislopeslinp4.append(0)
        pslopes5.append(0)
        pdeath5.append(0)
        scislopeslin5.append(0)
        scislopeslinr5.append(0)
        scislopeslinp5.append(0)
    elif len(indices2) == 2:
        pslopes1.append(b[0])
        pdeath1.append(TOD[1])
        maglinear1 = scipy.stats.linregress(xs[indices2[0]:indices2[1]-1], ys[indices2[0]:indices2[1]-1])
        scislopeslin1.append(maglinear1.slope)
        scislopeslinr1.append(maglinear1.rvalue**2)
        scislopeslinp1.append(maglinear1.pvalue)
        pslopes2.append(b[1])
        pdeath2.append(TOD[2])
        maglinear2 = scipy.stats.linregress(xs[indices2[1]:len(xs)], ys[indices2[1]:len(ys)])
        scislopeslin2.append(maglinear2.slope)
        scislopeslinr2.append(maglinear2.rvalue**2)
        scislopeslinp2.append(maglinear2.pvalue)
        pslopes3.append(0)
        pdeath3.append(0)
        scislopeslin3.append(0)
        scislopeslinr3.append(0)
        scislopeslinp3.append(0)
        pslopes4.append(0)
        pdeath4.append(0)
        scislopeslin4.append(0)
        scislopeslinr4.append(0)
        scislopeslinp4.append(0)
        pslopes5.append(0)
        pdeath5.append(0)
        scislopeslin5.append(0)
        scislopeslinr5.append(0)
        scislopeslinp5.append(0)
    elif len(indices2) == 3:
        pslopes1.append(b[0])
        pdeath1.append(TOD[1])
        maglinear1 = scipy.stats.linregress(xs[indices2[0]:indices2[1]-1], ys[indices2[0]:indices2[1]-1])
        scislopeslin1.append(maglinear1.slope)
        scislopeslinr1.append(maglinear1.rvalue**2)
        scislopeslinp1.append(maglinear1.pvalue)
        pslopes2.append(b[1])
        pdeath2.append(TOD[2])
        maglinear2 = scipy.stats.linregress(xs[indices2[1]:indices2[2]-1], ys[indices2[1]:indices2[2]-1])
        scislopeslin2.append(maglinear2.slope)
        scislopeslinr2.append(maglinear2.rvalue**2)
        scislopeslinp2.append(maglinear2.pvalue)
        pslopes3.append(b[2])
        pdeath3.append(TOD[3])
        maglinear3 = scipy.stats.linregress(xs[indices2[2]:len(xs)], ys[indices2[2]:len(ys)])
        scislopeslin3.append(maglinear3.slope)
        scislopeslinr3.append(maglinear3.rvalue**2)
        scislopeslinp3.append(maglinear3.pvalue)
        pslopes4.append(0)
        pdeath4.append(0)
        scislopeslin4.append(0)
        scislopeslinr4.append(0)
        scislopeslinp4.append(0)
        pslopes5.append(0)
        pdeath5.append(0)
        scislopeslin5.append(0)
        scislopeslinr5.append(0)
        scislopeslinp5.append(0)
    elif len(indices2) == 4:
        pslopes1.append(b[0])
        pdeath1.append(TOD[1])
        maglinear1 = scipy.stats.linregress(xs[indices2[0]:indices2[1]-1], ys[indices2[0]:indices2[1]-1])
        scislopeslin1.append(maglinear1.slope)
        scislopeslinr1.append(maglinear1.rvalue**2)
        scislopeslinp1.append(maglinear1.pvalue)
        pslopes2.append(b[1])
        pdeath2.append(TOD[2])
        maglinear2 = scipy.stats.linregress(xs[indices2[1]:indices2[2]-1], ys[indices2[1]:indices2[2]-1])
        scislopeslin2.append(maglinear2.slope)
        scislopeslinr2.append(maglinear2.rvalue**2)
        scislopeslinp2.append(maglinear2.pvalue)
        pslopes3.append(b[2])
        pdeath3.append(TOD[3])
        maglinear3 = scipy.stats.linregress(xs[indices2[2]:indices2[3]-1], ys[indices2[2]:indices2[3]-1])
        scislopeslin3.append(maglinear3.slope)
        scislopeslinr3.append(maglinear3.rvalue**2)
        scislopeslinp3.append(maglinear3.pvalue)
        pslopes4.append(b[3])
        pdeath4.append(TOD[4])
        maglinear4 = scipy.stats.linregress(xs[indices2[3]:len(xs)], ys[indices2[3]:len(ys)])
        scislopeslin4.append(maglinear4.slope)
        scislopeslinr4.append(maglinear4.rvalue**2)
        scislopeslinp4.append(maglinear4.pvalue)
        pslopes5.append(0)
        pdeath5.append(0)
        scislopeslin5.append(0)
        scislopeslinr5.append(0)
        scislopeslinp5.append(0)
    elif len(indices2) == 5:
        pslopes1.append(b[0])
        pdeath1.append(TOD[1])
        maglinear1 = scipy.stats.linregress(xs[indices2[0]:indices2[1]-1], ys[indices2[0]:indices2[1]-1])
        scislopeslin1.append(maglinear1.slope)
        scislopeslinr1.append(maglinear1.rvalue**2)
        scislopeslinp1.append(maglinear1.pvalue)
        pslopes2.append(b[1])
        pdeath2.append(TOD[2])
        maglinear2 = scipy.stats.linregress(xs[indices2[1]:indices2[2]-1], ys[indices2[1]:indices2[2]-1])
        scislopeslin2.append(maglinear2.slope)
        scislopeslinr2.append(maglinear2.rvalue**2)
        scislopeslinp2.append(maglinear2.pvalue)
        pslopes3.append(b[2])
        pdeath3.append(TOD[3])
        maglinear3 = scipy.stats.linregress(xs[indices2[2]:indices2[3]-1], ys[indices2[2]:indices2[3]-1])
        scislopeslin3.append(maglinear3.slope)
        scislopeslinr3.append(maglinear3.rvalue**2)
        scislopeslinp3.append(maglinear3.pvalue)
        pslopes4.append(b[3])
        pdeath4.append(TOD[4])
        maglinear4 = scipy.stats.linregress(xs[indices2[3]:indices2[4]-1], ys[indices2[3]:indices2[4]-1])
        scislopeslin4.append(maglinear4.slope)
        scislopeslinr4.append(maglinear4.rvalue**2)
        scislopeslinp4.append(maglinear4.pvalue)
        pslopes5.append(b[4])
        pdeath5.append(TOD[5])
        maglinear5 = scipy.stats.linregress(xs[indices2[4]:len(xs)], ys[indices2[4]:len(ys)])
        scislopeslin5.append(maglinear5.slope)
        scislopeslinr5.append(maglinear5.rvalue**2)
        scislopeslinp5.append(maglinear5.pvalue)
        

data1 = np.column_stack((freqnew, mznew[peak_idx], xdftnew[peak_idx], pslopes1, pslopes2, pslopes3, pslopes4, pslopes5, pdeath1, pdeath2, pdeath3, pdeath4, pdeath5, scislopeslin1, scislopeslin2, scislopeslin3, scislopeslin4, scislopeslin5, scislopeslinr1, scislopeslinr2, scislopeslinr3, scislopeslinr4, scislopeslinr5, scislopeslinp1, scislopeslinp2, scislopeslinp3, scislopeslinp4, scislopeslinp5))
slopesdf1 = pd.DataFrame(data1, columns = ['Frequency', 'm/z', 'Intensity', 'P Slope1', 'P Slope2', 'P Slope3', 'P Slope4', 'P Slope5', 'P Death1', 'P Death2', 'P Death3', 'P Death4', 'P Death5', 'SciLin Slope1', 'SciLin Slope2', 'SciLin Slope3', 'SciLin Slope4', 'SciLin Slope5', 'SciLin R-squared1', 'SciLin R-squared2', 'SciLin R-squared3', 'SciLin R-squared4', 'SciLin R-squared5', 'SciLin P-value1', 'SciLin P-value2', 'SciLin P-value3', 'SciLin P-value4', 'SciLin P-value5']) 
slopesdf1 = slopesdf1.sort_values('m/z')
slopesdf1['Scan'] = os.path.splitext(os.path.basename(pkl_files))[0]
slopesdf1['Ions'] = np.size(mznew[peak_idx])
slopesdf1.to_pickle(pkl_files+"_stori.pkl", compression = 'xz')


finish = tm.perf_counter()

print(f"Completed in {finish-start} seconds")
# <img src="https://github.com/clelandtp/Theropod/blob/39a582723d670a390fbf363506a01e0279c015b1/logo.png" width="500"/> 
# Theropod
Theropod STORI analysis tools

### Data acquisition with external digital oscilloscope
Either the pico-repeat-block-triggeron3V-queuing.py or pico-repeat-block-triggeron3V-queuing.ipynb can be used to collect transients. Updates to the storage folder name is required and must be unique. To use the scripts, the PicoSDK must be installed.

### Transient processing
Two versions of the the transient processor are available for single file processing (e.g., for job array on clusters) or a multiprocessor version that will process all *.ftr files in a directory.

#### For single file processing
```sh
$ python theropod_slope_single_file.py <FILENAME.ftr>
```

#### For multi file processing
```sh
$ python theropod_slope_multi_file.py
```

### Theropod Calibration factor calculation
Use the theropod_slope_calib_filter.py script to generate calibration files: one myoglobin file labeled with "myo" and one carbonic anhydrase 2 labeled with "cah2" are required in a single folder.
```sh
$ python theropod_slope_calib_filter.py <FILENAME>
```
In RStudio, theropod_calibration.R can be directly run after setting the working directory containing the myo and cah2 files. A KDE plot with the conversion factor will be output.

### Theropod SlopeComb Processing
To filter and perform charge assignment for the pkl files theropod.py is used. Command line parameters can be set as follows:
  
  --charge CHARGE       Set the number of charges to evaluate; Default = 60
  
  --neighbor NEIGHBOR   Set the number of either charge or isotope neighbors; Default = 10
  
  --tod TOD             Set the minimum ion life time; Default = 300
  
  --rsquare RSQUARE     Set the r-squared value; Default = 0.99
  
  --slopecal SLOPECAL   Set the KDE slope calibration value; Default = 2.292161
  
  --minmass MINMASS     Set the minimum mass for the final neutral mass spectrum; Default = 5000
  
  --maxmass MAXMASS     Set the maximum mass for the final neutral mass spectrum; Default = 75000
  
  --chargecorrecttype CHARGECORRECTTYPE     Set the charge correction type 1 = isotopes, 2 = charge neighbors; Default = 1
  
  --spectrumtype SPECTRUMTYPE     Set the spectrum type 1 = kernal density estimation, 2 = mass histogram; Default = 1

Tables containing the filtered data and final mass calculation will be exported. Additional plots of SlopeComb vs m/z; SlopeComb vs Mass; Slope Choice; m/z vs Final Charge; Final Charge Histogram; a representative m/z mass spectrum from 600-1200 m/z; ion count histogram. Additionally an interactive neutral mass spectrum is output in html format. The neutral mass spectrum and a M+H spectrum are exported as csv and mzML files.

The location of mzR installed through RStudio must be updated for the version of R used (line 4 of theropod_mzML.R) 

```sh
$ python theropod.py <FILENAME> --charge CHARGE --neighbor NEIGHBOR --tod TOD --rsquare RSQUARE --slopecal SLOPECAL --minmass MINMASS --maxmass MAXMASS --chargecorrecttype CHARGECORRECTTYPE --spectrumtype SPECTRUMTYPE
```
Minimum commandline
```sh
$ python theropod.py <FILENAME>
```
### Example STORI data
The example STORI pkl files are from carbonic anhydrase 2.

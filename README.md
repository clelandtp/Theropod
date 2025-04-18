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
To filter and perform charge assignment for the pkl files theropod_v1_0_0.py is used. Prior to full processing the calibration factor must be updated (Line 156). Settings for td (line 68), r-squared (line 69), charge shift (line 161) can be modified. Default settings for tod is 300 ms and r-squared is 0.99. Charge shift should be adjusted depending on expected multi-ion events and protein size. Tables containing the filtered data and final mass calculation will be exported. Additional plots of SlopeComb vs m/z; SlopeComb vs Mass; Slope Choice; m/z vs Final Charge; Final Charge Histogram; a representative m/z mass spectrum from 600-1200 m/z. Additionally an interactive neutral mass spectrum is output in html format. The neutral mass spectrum and a M+H spectrum are exported as csv and mzML files.

The location of mzR installed through RStudio must be updated for the version of R used (line 4 of theropod_mzML.R) 

```sh
$ python theropod_v1_0_0.py <FILENAME>
```

### Example STORI data
The example STORI pkl files are from carbonic anhydrase 2.

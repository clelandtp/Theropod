#
# Copyright (C) 2018 Pico Technology Ltd. See LICENSE file for terms.
#Copyright © 2019 Pico Technology Ltd.
#Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.
#THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#Modified for Theropod data acquistion by Timothy P. Cleland 2024

import ctypes
import numpy as np
import pandas as pd
from picosdk.ps4000a import ps4000a as ps
from picosdk.functions import adc2mV, assert_pico_ok
import time
import os
from queue import Queue
from threading import *

# Create chandle and status ready for use
chandle = ctypes.c_int16()
status = {}

# Open 4000 series PicoScope
# Returns handle to chandle for use in future API functions
status["openunit"] = ps.ps4000aOpenUnit(ctypes.byref(chandle), None)

try:
    assert_pico_ok(status["openunit"])
except: # PicoNotOkError:
    powerStatus = status["openunit"]

    if powerStatus == 286:
        status["changePowerSource"] = ps.ps4000aChangePowerSource(chandle, powerStatus)
    else:
        raise

    assert_pico_ok(status["changePowerSource"])

enabled = 1
disabled = 0
analogue_offset = 0.0
    
    
# Set up channel A
# handle = chandle
# channel = PS4000A_CHANNEL_A = 0
# enabled = 1
# coupling type = PS4000A_DC = 1
# range = PS4000A_5V = 8
# analogue offset = 0 V
channel_range = 8
status["setChA"] = ps.ps4000aSetChannel(chandle,
                                        ps.PS4000A_CHANNEL['PS4000A_CHANNEL_A'],
                                        enabled,
                                        ps.PS4000A_COUPLING['PS4000A_DC'],
                                        channel_range,
                                        analogue_offset)
assert_pico_ok(status["setChA"])

#set directory name for output
os.mkdir('20250115_pierceintactprotein_3ngul_0-3ulmin_0-5MIT_500ns_10000scans_240k2')

#~30 min @ 240k 2250 buffers
#~60 min @ 240k 4500 buffers
#~120 min @ 120k 9000 buffers
numBuffersToCapture = 10000

#240k 768 ms
#15k 48 ms
#60k 192 ms
#120k 384 ms
transient = 768
transsamp = np.round(transient/0.0005).astype(int)

#print(delay)
print(transsamp)

que = Queue()
def acquire(i):
    for i in range(1, i+1):          # Write each DataFrame to separate ftr file
        # Set up single trigger
        # handle = chandle
        # enabled = 1
        # source = PS4000a_CHANNEL_A = 0
        # threshold = 19660 ADC counts
        # direction = PS4000a_RISING = 2
        # delay = 57164 full scan samples with 0.5 ms MIT to only waveform
        #delay = 56683 full scan samples with 0.5 MIT to OT trigger
        #delay = 86748 HCD with 0.5 MIT to only waveform
        #delay = 63404 full scan with 5.0 MIT to only waveform
        #delay = 483 for 725 ns sampling
        #delay = 700 for 500 ns sampling
        # auto Trigger = 5000 ms
        status["trigger"] = ps.ps4000aSetSimpleTrigger(chandle, 1, 0, 19660, 2, 700, 0)
        assert_pico_ok(status["trigger"])

        # Set number of pre and post trigger samples to be collected
        preTriggerSamples = 0
        postTriggerSamples = transsamp
        maxSamples = preTriggerSamples + postTriggerSamples

        # Get timebase information
        # handle = chandle
        # timebase = 57 = timebase = 725 ns
        # timebase = 39 = timebase = 500 ns
        # noSamples = maxSamples
        # pointer to timeIntervalNanoseconds = ctypes.byref(timeIntervalns)
        # pointer to maxSamples = ctypes.byref(returnedMaxSamples)
        # segment index = 0
        timebase = 39
        timeIntervalns = ctypes.c_float()
        returnedMaxSamples = ctypes.c_int32()
        oversample = ctypes.c_int16(1)
        status["getTimebase2"] = ps.ps4000aGetTimebase2(chandle, timebase, maxSamples, ctypes.byref(timeIntervalns), ctypes.byref(returnedMaxSamples), 0)
        assert_pico_ok(status["getTimebase2"])
       
        # Run block capture
        # handle = chandle
        # number of pre-trigger samples = preTriggerSamples
        # number of post-trigger samples = PostTriggerSamples
        # time indisposed ms = None (not needed in the example)
        # segment index = 0
        # lpReady = None (using ps4224aIsReady rather than ps4224aBlockReady)
        # pParameter = None
        status["runBlock"] = ps.ps4000aRunBlock(chandle, preTriggerSamples, postTriggerSamples, timebase, None, 0, None, None)
        assert_pico_ok(status["runBlock"])

        # Check for data collection to finish using ps4224a
        ready = ctypes.c_int16(0)
        check = ctypes.c_int16(0)
        while ready.value == check.value:
            status["isReady"] = ps.ps4000aIsReady(chandle, ctypes.byref(ready))

        # Create buffers ready for assigning pointers for data collection
        bufferAMax = (ctypes.c_int16 * maxSamples)()
        bufferAMin = (ctypes.c_int16 * maxSamples)()

        # Set data buffer location for data collection from channel A
        # handle = chandle
        # source = PS4000a_CHANNEL_A = 0
        # pointer to buffer max = ctypes.byref(bufferAMax)
        # pointer to buffer min = ctypes.byref(bufferAMin)
        # buffer length = maxSamples
        # segementIndex = 0
        # mode = PS4000A_RATIO_MODE_NONE = 0
        status["setDataBuffersA"] = ps.ps4000aSetDataBuffers(chandle, 0, ctypes.byref(bufferAMax), ctypes.byref(bufferAMin), maxSamples, 0 , 0)
        assert_pico_ok(status["setDataBuffersA"])

        # create overflow loaction
        overflow = ctypes.c_int16()
        # create converted type maxSamples
        cmaxSamples = ctypes.c_int32(maxSamples)

        # Retried data from scope to buffers assigned above
        # handle = chandle
        # start index = 0
        # pointer to number of samples = ctypes.byref(cmaxSamples)
        # downsample ratio = 0
        # downsample ratio mode = PS4000a_RATIO_MODE_NONE
        # pointer to overflow = ctypes.byref(overflow))
        status["getValues"] = ps.ps4000aGetValues(chandle, 0, ctypes.byref(cmaxSamples), 0, 0, 0, ctypes.byref(overflow))
        assert_pico_ok(status["getValues"])
        que.put(bufferAMax)

#Change foldername 
def writer(j):
    for j in range(1,j+1):
        start = time.perf_counter()
        cols=['Channel A']
        picospec = que.get()
        spectrum = pd.DataFrame(np.transpose(picospec), columns = cols)
        spectrum.to_feather(r'20250115_pierceintactprotein_3ngul_0-3ulmin_0-5MIT_500ns_10000scans_240k2\scan-'+ str(j) +'.ftr', compression = 'zstd')
        print('scan-'+str(j))
        finish = time.perf_counter()
        print(f"Completed in {finish-start} seconds")
        

t1=Thread(target=acquire, args=(numBuffersToCapture,))
t2=Thread(target=writer, args=(numBuffersToCapture,))
t1.start()
t2.start()

   
t1.join()
t2.join()

# Stop the scope
# handle = chandle
status["stop"] = ps.ps4000aStop(chandle)
assert_pico_ok(status["stop"])

# Close unitDisconnect the scope
# handle = chandle
status["close"] = ps.ps4000aCloseUnit(chandle)
assert_pico_ok(status["close"])

# display status returns
print(status)


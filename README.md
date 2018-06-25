
# filtershift.cc

Filtershift.cc is a c++ module for fMRI preprocessing developed by QNL in the Cognitive Neuroscience Division of the Columbia University Medical Center Neurology Department.  Filtershift uses the principles of digital signal processing and ideal signal reconstruction to optimally resample a time series at a constant offset.  It also (optionally) filters the data using an optimized kaiser-window LPF.

## Reference:

[Parker, David, Xueqing Liu, and Qolamreza R. Razlighi. "Optimal slice timing correction and its interaction with fMRI parameters and artifacts." Medical Image Analysis 35 (2017): 434-445.](http://dx.doi.org/10.1016/j.media.2016.08.006)

## Help

### Application

Filtershift is a slice timing correction routine.  If you use Filtershift in your preprocessing pipeline, you should remove all other slice timing correction modules (also known as temporal realignment).    

## Notes

### "Slice Timing" file VS "Slice Order" File

#### Slice Order:
This is the typical slice order format. "Row" indicates the order in which the slices are acquired (Slice acquired 1st, slice acquired 2nd, etc), "val" indicates the physical slice index that was acquired.

  Row:  Val:
  1     1
  2     3
  3     5
  4     7  
  5     2
  6     4
  7     6
  8     8

 In this example, row 2 - val 3 means slice 3 was acquired 2nd.  row 5 - val 2 means slice 2 was acquired 5th, and so on

#### Slice Timing
This is the typical slice timing format.  "Row" indicates the physical slice index, and "Val" indicates the time at which that slice was acquired. For example, row 1 is the time of slice 1's acquisition, row 2 is the time of slice 2's acquisition, etc

  Row:  Val:
  1     0.000
  2     0.500
  3     0.033
  4     0.533  
  5     0.066
  6     0.566
  7     0.100
  8     0.600

Now in this case, row 2 val 0.500 means slice 2 was acquired 0.500 seconds after the first slice acquired in the volume.  Row 5 val 0.066 means slice 5 was acquired 0.066 seconds after the start of the first slice slice acquired in the volume.

### Multiband

#### Slice Timing
Multiband works fine with a slice timing file.  You just manually set the time that each slice was acquired in a way that accounts for the simultaneous acquisition of multiple slices.  For example:

  Row:  Val:
  1     0.000
  2     0.666
  3     0.333
  4     1.000  
  5     0.000
  6     0.666
  7     0.333
  8     1.000
  
  This will correctly account for slices 1, 4, and 7 being acquired simultaneously at time 0.000s, slices 3 and 6 being acquired simultaneously at time 0.500s, and slices 2, 5, and 8 being acquired simultaneously at time 1.000s.
  
#### Slice Order
The program can detect a multiband slice order file.  The format of the multiband slice order file is fundamentally different from the slice order file for a sequential volume.  Please read the descriptions of each, and ensure that you fully understand the differences.

 Multiband slice order files use the following format, which means slice "row" was acquired "val"

  Row:  Val:
  1     1
  2     3
  3     2
  4     4  
  5     1
  6     3
  7     2
  8     4

Now in this case, row 2 val 3 means slice 2 was acquired 3rd.  Row 5 val 1 means slice 5 was acquired 1st.  Note that this example is equivalent to the example given for the multiband slice Timing file.  You would basically put the order in if you were too lazy to just calculate the time shifts yourself.  Yes, I'm looking at you. 

### Basic Installation

#### Binary File
You can simply download the binary executable file "filtershift" and add it to a directory of your choice.  Then, modify your .bashrc file so your PATH variable includes that directory.  

#### Manual Compile
You can manually compile the source code yourself, however this requires that you have fsl 5.0.7 or later installed on your machine.
##### 1. Navigate to the "src" folder in your FSL directory ('/usr/local/fsl/5.0.7' or similar)
##### 2. In the "src" folder, create a new directory called "filtershift"
##### 3. Download the associated \*.h and \*.cc files to this directory
##### 4. in the FSL 5.0.7 directory, call the "build" file in terminal, with the argument "filtershift":
```shell
./build filtershift
```
Note: you will need superuser access to do this successfully.
##### 5. you should now be able to run ./filtershift from a new terminal window, providing your FSL paths are set up correctly.


## Using filtershift
```shell
filtershift --in <InputFile> --tr <TR> [options]

Compulsory arguments (You MUST set one or more of):
	--in		filename the input image to perform STC on

	--TR		Set the TR of the original fMRI data in seconds


filtershift --in <InputFile> --tr <TR> [options]
	 If no other options are set, this assumes ascending
	 Slice order with no interleave.

	 To use with a slice Order file, the TR must be
	 specified, and slices will be shifted according to
	 a fractional amount of the TR.  By default, we align
	 to the first slice acquired in the TR.

	 To use a slice timing file, no TR is required,
	 simply provide the desired shift in seconds as a
	 Column vector, saved in a text file.  Slice
	 Shifting is independent of all other slices,
	 so multiple slices can be shifted the same amount
	 (To correct multi-band acquisition, for example)


Compulsory arguments (You MUST set one or more of):
	-i,--in	 filename the input image to perform STC on


Optional arguments (You may optionally specify one or more of):
	-h,--help	 display this message

	--TR	 Set the TR of the original fMRI data in seconds

	--itl	 set the interleave parameter, or how many slices are
			 incremented between acquisitions
			 1 = sequential acquisition
			 2 = even/odd acquisition (acquire every second slice)
			 3 = acquire every third slice, etc
			 Leaving this blank will assume bottom up sequential
			 acquisition

	-o,--out	 Specify an output file name - all working directories
			 will be created in the parent directory specified
			 here. Leave blank to run in the parent directory
			 of <InputFile>

	-s,--start	 Set the starting slice - The slice that was acquired
			 first in the sequence. Default is slice 1, the bottom
			 most slice. This starts the interleave from that slice.
			 If your interleave parameter is '1' and your starting
			 slice is '3', your slice acquisition sequence will be
			 modeled as:
				 3
				 5
				 7
				 9...

	-d,--direction	 value 1 or -1.  Set the direction of slice 
			 acquisition.
			 1: ascending slice acquisition:(1,3,5,7,9...)
			-1: descending slice acquisition: (9,7,5,3,...)

	--order		 Slice Order File.  This file is the order in which
			 each slice was acquired. each row represents the
			 order in which that slice was acquired. For example,
			 '1' in the first row means that slice 1 was acquired
			 first. '20' in the second row means that slice 20 was
			 acquired 2nd. If present, all interleave parameters
			 are ignored, and slices are shifted using the slice
			 order file. we refer to the bottom slice in the image
			 as slice 1, not slice 0

	--cf		 Set the cutoff frequency of the lowpass filter in Hz.
			 Note that by default, this is set to 0.21 Hz.

	--timing	 Slice Timing File.  This file is the time at which
			 each slice was acquired relative to the first slice.
			 each row represents the time at which that slice was
			 acquired. For example, '0' in the first row means
			 that slice 1 was acquired first, and will be shifted 0
			 seconds. '0.5' in the second row means that slice 2
			 was acquired 0.5 seconds after the first slice, and
			 will be shifted 0.5s. If present, all interleave
			 parameters are ignored, and slices are shifted using
			 the slice timing file. If you put both a slice timing
			 and a slice order, the program will yell at you and
			 refuse to run.

	-r,--reference	 Set the Reference slice
			 This is the slice the data is aligned to.
			 Default is the first slice

	--lpf	 Only Run the Lowpass Filter, do not
			 Preform slice timing correction
```
## About QNL

[QNL Homepage](http://www.columbia.edu/cu/qnl/)



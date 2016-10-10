#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>


#include "miscmaths/optimise.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "utils/options.h"
#include "miscmaths/kernel.h"
#include "Window.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace Utilities;

string title="\n\nFilterShift \nQNL's tool for optimal slice timing correction\n\
#WeRememberHarambe\n";
string examples="filtershift --in <InputFile> --tr <TR> [options]\
\n\t If no other options are set, this assumes ascending\
\n\t Slice order with no interleave.\n\
\n\t To use with a slice Order file, the TR must be\
\n\t specified, and slices will be shifted according to\
\n\t a fractional amount of the TR.  By default, we align\
\n\t to the first slice acquired in the TR.\n\
\n\t To use a slice timing file, no TR is required,\
\n\t simply provide the desired shift in seconds as a\
\n\t Column vector, saved in a text file.  Slice\
\n\t Shifting is independent of all other slices,\
\n\t so multiple slices can be shifted the same amount\
\n\t (To correct multi-band acquisition, for example)\n";

Option<string> input(string("-i,--in"), string(""),string(
" filename the input image to perform STC on\n"),
			  true, requires_argument);

Option<float> TR(string("--TR"), 2.0,string(
" Set the TR of the original fMRI data in seconds\n"),
			  false, requires_argument);

Option<int> interleave(string("--itl"), 0, string(
" set the interleave parameter, or how many slices are\
\n\t\t\t incremented between acquisitions\
\n\t\t\t 1 = sequential acquisition\
\n\t\t\t 2 = even/odd acquisition (acquire every second slice)\
\n\t\t\t 3 = acquire every third slice, etc\
\n\t\t\t Leaving this blank will assume bottom up sequential\
\n\t\t\t acquisition\n"), 
			 false, requires_argument);

Option<string> out(string("-o,--out"), string(""),string(
	   " Specify an output file name - all working directories\
\n\t\t\t will be created in the parent directory specified\
\n\t\t\t here. Leave blank to run in the parent directory\
\n\t\t\t of <InputFile>\n"), 
			 false, requires_argument);

Option<int> start(string("-s,--start"), 0,string(
	   " Set the starting slice - The slice that was acquired\
\n\t\t\t first in the sequence. Default is slice 1, the bottom\
\n\t\t\t most slice. This starts the interleave from that slice.\
\n\t\t\t If your interleave parameter is '1' and your starting\
\n\t\t\t slice is '3', your slice acquisition sequence will be\
\n\t\t\t modeled as:\
\n\t\t\t\t 3\
\n\t\t\t\t 5\
\n\t\t\t\t 7\
\n\t\t\t\t 9...\n"),
			  false, requires_argument);
			  
Option<int> direction(string("-d,--direction"), 1,string(
		  " value 1 or -1.  Set the direction of slice \
\n\t\t\t acquisition.\
\n\t\t\t 1: ascending slice acquisition:(1,3,5,7,9...)\
\n\t\t\t-1: descending slice acquisition: (9,7,5,3,...)\n"), 
			 false, requires_argument);

Option<string> order(string("--order"), string(""),string(
     "\t Slice Order File.  This file is the order in which\
\n\t\t\t each slice was acquired. each row represents the\
\n\t\t\t order in which that slice was acquired. For example,\
\n\t\t\t '1' in the first row means that slice 1 was acquired\
\n\t\t\t first. '20' in the second row means that slice 20 was\
\n\t\t\t acquired 2nd. If present, all interlave parameters\
\n\t\t\t are ignored, and slices are shifted using the slice\
\n\t\t\t order file. we refer to the bottom slice in the image\
\n\t\t\t as slice 1, not slice 0\n"),
			 false, requires_argument);

Option<float> cf(string("--cf"),0.21,string(
	 "\t Set the cutoff frequency of the lowpass filter in Hz.\
\n\t\t\t Note that by default, this is set to 0.21 Hz.\n"),
				  false,requires_argument);

Option<string> timing(string("--timing"), string(""),string(
	   " Slice Timing File.  This file is the time at which\
\n\t\t\t each slice\ was acquired relative to the first slice.\
\n\t\t\t each row represents the time at which that slice was\
\n\t\t\t acquired. For example, '0' in the first row means\
\n\t\t\t that slice 1 was acquired first, and will be shifted 0\
\n\t\t\t seconds. '0.5' in the second row means that slice 2\
\n\t\t\t was acquired 0.5 seconds after the first slice, and\
\n\t\t\t will be shifted 0.5s. If present, all interleave\
\n\t\t\t parameters are ignored, and slices are shifted using\
\n\t\t\t the slice timing file. If you put both a slice timing\
\n\t\t\t and a slice order, the program will yell at you and\
\n\t\t\t refuse to run.\n"),
			 false, requires_argument);

Option<int> ref(string("-r,--reference"), 1,string(
	   " Set the Reference slice\
\n\t\t\t This is the slice the data is aligned to.\
\n\t\t\t Default is the first slice\n"),
			 false, requires_argument);

Option<int> lpf(string("--lpf"), 1,string(
	   " Only Run the Lowpass Filter, do not\
\n\t\t\t Preform slice timing correction\n"),
			 false, no_argument);
			      
  Option<bool> help(string("-h,--help"), false,
		  string(" display this message\n"),
		  false, no_argument);


inline bool exists_test3 (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


void filter_timeseries(ColumnVector *timeseries, std::vector<float> *FIR, int shift,int skip,int slice)
{
	
// 		9/27/16 - modified Kaiser Window resampling algorithm and convolution filtering routine.
// 		Now matches output from old code almost perfectly (10e-3 error)
	

	int lenT = timeseries->Nrows();
	int lenF = FIR->size();
	
	if ( lenF/skip >= lenT )
	{
		std::cout<<"Filter Order too high.  There aren't enough time points in your image."<< std::endl;
		return;
	}
		
	std::vector<int> SamplePoints;
	
	for (int i=1;i<=lenF-skip;i+=skip)
	{
		SamplePoints.insert(SamplePoints.end(),i);
	}
	
	int firLen=SamplePoints.size();	
	std::vector<float> pFIR;
	pFIR.assign(FIR->begin(),FIR->end());
	std::vector<float> padd;
	padd.reserve(std::abs(shift));	
	int ModSample=0;
	
	// If the shift if positive (shifting the signal to the right), then we want to DELAY the filter, add zeros to the END (Right hand side)		
	if (shift>0)
	{

		padd.assign(std::abs(shift),pFIR.back());
		pFIR.insert (pFIR.end(),padd.begin(),padd.end());
		ModSample=shift;
	}
	
	else if (shift<0)
	{
		padd.assign(std::abs(shift),pFIR.front());

		padd.insert(padd.end(),pFIR.begin(),pFIR.end());
		pFIR=padd;
		
		ModSample=0;
		
	}	
	
	ColumnVector FIR_down_shift;
	ColumnVector FIR_down;
	FIR_down_shift.ReSize(firLen);
	FIR_down.ReSize(firLen);
	lenF=FIR_down_shift.Nrows();	
	firLen=1;
	
	for (int i = 0; i< SamplePoints.size(); i++)
	{
		FIR_down(firLen)=FIR->operator[](SamplePoints[i]);
		SamplePoints[i]+=ModSample;
		FIR_down_shift(firLen)=pFIR[SamplePoints[i]];
		firLen+=1;
	}

	ColumnVector filtered;
	filtered.ReSize(lenT);	
	int startT = floor(lenF/2);	
	int maxT = lenT-lenF-1;
	float FiltSum=0;

	for (int i = 0; i<maxT; i++)
	{
		FiltSum=0;

		for (int f = 1; f<=lenF; f++)
		{
			FiltSum+=FIR_down_shift(f)*timeseries->operator()(i+f);
		}
		
		filtered(i+startT)=FiltSum;		
	}

	ColumnVector filtered2;
	filtered=filtered.Reverse();
	filtered2=filtered;
	
	for (int i = 0; i<maxT; i++)
	{
		FiltSum=0;
		
		for (int f = 1; f<=lenF; f++)
		{
			FiltSum+=FIR_down(f)*filtered(i+f);
		}
		
		filtered2(i+startT)=FiltSum;
	}
	
	filtered=filtered2.Reverse();
	*timeseries=filtered;
}

void make_timings(Matrix *timings, Matrix *orders, int zs)
{
	if ( timing.set() )
	{
		
		if ( ref.set() )
		{
			std::cout<<"When using a Slice Timing file, the times in the file supercede all other settings.  The reference slice specified will be ignored"<<std::endl;
			std::cout<<"If you wish to alighn data to a specific slice, please make that adjustment in the Slice Timing File, or omit the slice timing file."<<std::endl;
			std::cout<<"Or use a Slice ORDER file, and specify a reference slice that way.  There's one clear option here that involves the least amount of work."<<std::endl;
		}
		// Slice Timing File, deftault Reference Slice, Tested 9/22/16 - Shifting Success, Slice Order Not
		// Slice Timing File, Custom Reference Slice, Tested 9/22/16 - Shifting Success, Slice Order Not.
		// Slice Order is Unnecessary, removing 9/23/16
		
		// Need To Test With Multiband Images
		
		Matrix TempTimings;
		TempTimings.ReSize(zs,1);
		
		try
		{
			*timings = read_ascii_matrix(timing.value(), zs, 1);
		}
		catch (...)
		{
			std::cout<<"Error Loading file "<<timing.value()<<std::endl;
			return;
		}
		
		
		TempTimings = read_ascii_matrix(timing.value(), zs, 1);
		if (timings->Nrows()!=zs)
		{
			std::cout<<"Slice timing file does not have the correct number of slices"<<std::endl;
			return;
		}
		
		// 9/23/16 - removed code that calculated slice order - it's wrong and unnecessary
		

		
	}	
	else if ( order.set() )
	{
		//Slice Order File, Default Reference Slice Tested 9/23/16 - Passed
		//Slice Order File, Custom Reference Slice Tested 9/23 - Passed		
		// 9/23/16 - Will not Run With Multiband, only slice timing file (User Burden, Deal With It)
		//Slice Order File, Default Reference Slice Tested 10/10/16 - Passed
		//Slice Order File, Custom Reference Slice Tested 10/10/16 - Passed		
		
		int tmx;
		float shift;
		
		try
		{
			*orders = read_ascii_matrix(order.value(), zs, 1);
		}
		catch (...)
		{
			std::cout<<"Error Loading file "<<order.value()<<std::endl;
			return;
		}
		tmx=orders->Maximum();
		if (orders->Nrows()!=zs)
		{
			std::cout<<"Slice order file does not have the correct number of slices"<<std::endl;
			return;
		}
		Matrix TimeList;
		TimeList.ReSize(zs,1);		
		shift=TR.value()*1.0/tmx;

		for ( int i=1; i<=zs; i++ )
		{
			TimeList(i,1)=(float) (i-1)*shift;
		}

		// 9/23/16 - Made This loop more efficient
		
		for (int j=1;j<=zs;j++)
		{			
			timings->operator()(orders->operator()(j,1),1)=TimeList(j,1);					
		}
		
		timings->operator-=(timings->operator()(ref.value(),1));
		
	}
	else
	{
		// Create Timing File Tested 9/22/16 - Succesful
		// With Reference Tested 9/22/16 - Succesful
		
		float dt;
		Matrix IntSeq;
		Matrix TimeList;
		IntSeq.ReSize(zs,1);		
		TimeList.ReSize(zs,1);
		dt=TR.value()/zs;
		int counter;
		counter=1;
		
		
		for ( int i=0; abs(i)<interleave.value(); i=i+1*direction.value() )
		{
			for ( int j =0; j<=floor((float) zs/interleave.value()); j++ ) 
			{
				if ((i)+(j*interleave.value())+1<=zs)
				{
					orders->operator()(counter,1)=(i)+(j*interleave.value())+1;
					TimeList((i)+(j*interleave.value())+1,1)=counter;
					counter++;
				}
			}
		}
		
		TimeList*=dt;
		TimeList-=TimeList(ref.value(),1);
		*timings=TimeList;
		

	}
}


int shift_volume()
{
	
	// 10-10-16 Added check for some combination of optional arguments
	// Now you idiots can't mess up your data by leaving out these settings.
	// ...but I'm sure you'll still manage to find a way >:(
	
	if (!TR.set()&&!order.set()&&!timing.set()&&!lpf.set())
	{
		std::cout<<"Must Specify either a TR, a slice Order file,a slice timing file,\n or indicate that you only wish to preform filtering"<<std::endl;
		return -1;
	}
	
	
	volume4D<float> timeseries;
	Matrix timings;
	Matrix orders;
	std::vector<float> FIR;
	
  if (input.set())
  {
	
	if (true) { cout << "Reading input volume" << endl; }
	read_volume4D(timeseries,input.value());

  } else if (out.set()) {
	cerr << "Must specify an input volume (-i or --in) to generate corrected data." << endl;
	return -1;
  }
  	

	
	timings.ReSize(timeseries.zsize(),1);
	orders.ReSize(timeseries.zsize(),1);	
	make_timings(&timings,&orders,timeseries.zsize());	
	
	int no_volumes = timeseries.tsize();
	int xx = timeseries.xsize();
	int yy = timeseries.ysize();
	int zz = timeseries.zsize();
	
	float cutoff=cf.value();
	float samplingrate=(float) zz/TR.value();
	float stopgain=-20;
	double transwidth=.1;
	
	window::window kaiser(cutoff,samplingrate,stopgain,transwidth);
	FIR=kaiser.get_fir();
	
	ColumnVector voxeltimeseries = timeseries.voxelts(1,1,1);
	ColumnVector fliptimeseries = voxeltimeseries.Reverse();
	ColumnVector cattimeseries;
	
	fliptimeseries=fliptimeseries.Rows(2,no_volumes-1);
	
	int lents=fliptimeseries.Nrows();
	int rangelh;
	int rangerh;
	
	if ( lents % 2 == 0 )
	{
		rangelh=lents/2;
		rangerh=lents/2-1;		
	}
	
	else
	{
		rangelh=ceil(lents/2);
		rangerh=floor(lents/2);
	}
	
	int cutLeft = lents-rangelh+2;
	int cutRight = cutLeft+no_volumes-1;
	
	if (cutRight-cutLeft != no_volumes-1)
	{
		return 1;
	}

	float span;
	float mn;
	float mn2;
	float span2;
	
	// 10/19/16 - "lpf" option now specifies JUST filtering.  This part sets
	// the slice delay to be zero for all of them (No shifting.  NO SHIFTING!)
	
	if (lpf.set())
	{
		for (int tm=1;tm<=zz;tm++)
		{
			timings(tm,1)=0;
		}
	}
	
	for (int slice=1; slice<=zz; slice++)
	{		
		
		for (int x_pos = 0; x_pos < xx; x_pos++)
		{
			
			for (int y_pos = 0; y_pos < yy; y_pos++)
			{
				voxeltimeseries = timeseries.voxelts(x_pos,y_pos,slice-1);
				mn=voxeltimeseries.Sum()/voxeltimeseries.Nrows();
				voxeltimeseries-=(mn);
				span=voxeltimeseries.Maximum()-voxeltimeseries.Minimum();
				
				// 9/27/16 - Filter changes mean and span - added code to maintain mean and span.
				// Also checks to see if the span is zero.  If so, do not filter.
				
				if ( span!=0 )
				{
					fliptimeseries = voxeltimeseries.Reverse();
					fliptimeseries=fliptimeseries.Rows(2,no_volumes-1);				
					cattimeseries=fliptimeseries.Rows(rangelh,lents)&voxeltimeseries&fliptimeseries.Rows(1,rangerh);				
					filter_timeseries(&cattimeseries, &FIR, (int)floor(timings(slice,1)*(float)samplingrate+0.5),zz,slice);
					cattimeseries=cattimeseries.Rows(cutLeft,cutRight);
					mn2=cattimeseries.Sum()/cattimeseries.Nrows();
					cattimeseries-=(mn2);
					span2=cattimeseries.Maximum()-cattimeseries.Minimum();
					cattimeseries/=span2;
					cattimeseries*=span;
					cattimeseries+=mn;					
					timeseries.setvoxelts(cattimeseries,x_pos,y_pos,slice-1);
				}
	
			}
		}
	}
	
	write_volume4D(timeseries,out.value());
  return 0;
}





int main (int argc,char** argv)
{
  
  OptionParser options(title, examples);

  try {
	options.add(help);
	options.add(input);
	options.add(TR);
	options.add(interleave);
	options.add(out);
	options.add(start);
	options.add(direction);
	options.add(order);
	options.add(cf);
	options.add(timing);
	options.add(ref);
	options.add(lpf);


	options.parse_command_line(argc, argv);

	if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
	  {
	options.usage();
	exit(EXIT_FAILURE);
	  }
		
	if ( input.unset()) 
	  {
	options.usage();
	cerr << endl 
		 << "--in or -i MUST be used." 
		 << endl;
	exit(EXIT_FAILURE);
	  }
	
  }  catch(X_OptionError& e) {
	options.usage();
	cerr << endl << e.what() << endl;
	exit(EXIT_FAILURE);
  } catch(std::exception &e) {
	cerr << e.what() << endl;
  } 

  if (!out.set())
  {
	string InputFile=input.value();
	string filename=InputFile.substr(InputFile.find_last_of( '/' ) + 1 );
	string::size_type n = filename.find_last_of('.');
	string directory=InputFile.substr(0,InputFile.find_last_of('/')+1);

	while (n!=string::npos)
	{
		filename=filename.substr(0,n);
		n=filename.find_last_of('.');
	}
	  out.set_value(directory+filename+ "_st.nii.gz");			
  }
  
  if (lpf.value()==0)
  {
	float cutoff=(float) 1.0/(2.0*TR.value())*0.8;
	std::ostringstream ss;
	ss << cutoff;
	std::string s(ss.str());
	lpf.set_value(s);
  }
  
  int retval = shift_volume();
  
  if (retval!=0) {
	cerr << endl << endl << "Error detected: try -h for help" << endl;
  }
	
	
  return retval;
}
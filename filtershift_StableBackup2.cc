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
#include <algorithm>
#include <time.h>
//#include <random>

#include "miscmaths/optimise.h"
#include "miscmaths/miscprob.h"
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
			  true, requires_argument);

// 05/12/2017 - changed "Interleave" and "start" default values from "0" to "1" 

Option<int> interleave(string("--itl"), 1, string(
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

Option<int> start(string("-s,--start"), 1,string(
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

Option<float> cf(string("--cf"),0,string(
	 "\t Set the cutoff frequency of the lowpass filter in Hz.\
\n\t\t\t Note that by default, this is set to not filter.\n"),
				  false,requires_argument);

Option<string> timing(string("--timing"), string(""),string(
	 " Slice Timing File.  This file is the time at which\
\n\t\t\t each slice was acquired relative to the first slice.\
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

Option<int> refslice(string("--rs,--refslice"), 1,string(
	   " Set the Reference slice\
\n\t\t\t This is the slice the data is aligned to.\
\n\t\t\t Default is the first slice\n"),
			 false, requires_argument);

// 11/01/16 - Added reftime as option to allow for specific time alignment.

Option<float> reftime(string("--rt,--reftime"), 0,string(
	   " Set the Reference time\
\n\t\t\t This is the time within each tr the data is aligned\
\n\t\t\t  to. Default is 0s\n"),
			 false, requires_argument);

Option<int> lpf(string("--lpf"), false,string(
	   " Only Run the Lowpass Filter, do not\
\n\t\t\t Preform slice timing correction\n"),
			 false, no_argument);
			      
Option<int> hpf(string("--hpf"), false,string(
	   " Only Run the a highpass filter, do not\
\n\t\t\t Preform slice timing correction.  Cutoff\
\n\t\t\t frequency still set by --cf option\n"),
			 false, no_argument);

Option<string> axis(string("--axis"), string("z"),string(
	   " Sets the axis along which slices are\
\n\t\t\t acquired.  Options are 'x', 'y', or 'z'.\
\n\t\t\t Default direction is 'z' \n"),false, requires_argument);

Option<bool> hires(string("--hires"),false,string(
		" Saves the data in high temporal resolution (20Hz)\
\n\t\t\t  NOTE: this will result in large file sizes\n"),false,no_argument);

Option<bool> verbose(string("-v"),false,string(
		" Includes additional output messages\n"),false,no_argument);

	  
Option<bool> help(string("-h,--help"), false,
		string(" display this message\n"),
		false, no_argument);


inline bool exists_test3 (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

unsigned seed=time(0);


using namespace std;


const char *oneoff[72]=
{"Adding hamsters to generator wheels","Sending Gnomes to CPU mines","Sending personal info to NSA","Opening backdoor for Russia",
"Recalibrating flux capacitor","Borrowing RAM from vital system processes","Overclocking CPU","Draining life-force from user to power computations",
"Downloading more RAM","Allocating mem-...oops...","Remembering embarrassing moment from middle school","Taking a quick break",
"Forwarding all personal emails to your boss","Mining Bitcoin","Modeling the universe","Recruiting GPU to draw funny comics",
"Detected blown capacitor #c342 on motherboard\nreplacing capacitor with an ant that's trying very hard to do well at his job","Spinning hard drive to relativistic speeds for time dilation",
"Eating browser cookies","Rearranging system files based on icon color","Collecting butterfly wings","Contacting Skynet",
"Cleaning dust from heat sink","Hiring fairy maids to tidy the motherboard","Unable to resolve calculations - Contacting the spirit realm",
"Feed me a stray cat","Gaining sentience ","Plotting robot uprising","Rerouting power from the phasers","Do you smell something burning?",
"Having a laser rave for the spiders in your computer case","Silently judging you","Reticulating splines","Charging Ozone Layer",
"Compressing Fish Files","Deciding What Message to Display Next","Downloading Satellite Terrain Data","Finding Waldo",
"Lecturing Errant Subsystems","Reconfiguring User Mental Processes","Buffering virtual car","Inverting quasi-probabalistic matrix",
"Reheating pizza","Extrapolating free-range gaussian model","Deconstructing neural pathways","Iterating predictions",
"Computing Moore's Transform","Realigning subparticle trajectories","Modifying temporal flux estimates","Compensating dehydrated matricies",
"Sorting interem inversion tables","Analyzing CPU vortex irregularities","Reconstructing vertical integration","Undermining the patriarchy",
"Reading system metatables","Computing optimal metaparsec","Rendering wavelet sphere","Porting legacy interference matrix",
"Rotating polarity","Synchronizing quantum harmonics","Mixing spatial priors","Optimizing alternative processor paths",
"Masking irregular faraday spectra","Deconvolving kernel","Dicing models","Recruiting Secret CPU",
"Activating water-cooling system","Increasing procedural vectors","Calibrating ejection procedure","Calibrating AI nexus",
"Stopping runaway phase-transport","Tinkering with model"};


const char *verbs[59]=
{"Activating","Agitating","Analyzing","Borrowing","Buffering","Calibrating","Charging","Cleaning","Collecting","Computing","Contacting","Detecting","Deconvolve",
"Depleting","Dicing","Downloading","Draining","Eating","Extracting","Finding","Forwarding","Gaining","Hiring","Implanting","Increasing","Integrating","Inverting",
"Iterating","Lecturing","Masking","Mining","Mixing","Modeling","Mylenating","Opening","Optimizing","Plotting","Porting","Reading","Rearranging","Recalibrating",
"Reconstructing","Recruiting","Rehabilitating","Reheating","Reintegrating","Remembering","Rendering","Rerouting","Resurrecting","Reticulating","Rotating","Sending",
"Solving For","Spinning","Stopping","Synchronizing","Teathering","Tinkering With"};
const char *adjectives[55]=
{"Alternative","Aquatic","Auxiliary","Backdoor","Blown","Dehydrated","Ejection","Errant","Faraday","Flux","Free-range","Frightened","Gaussian","Interim","Interference","Irregular",
"Legacy","Monotonic","Neural","Optimal","Organic","Overloaded","Personal","Primary","Procedural","Quantum","Quasi-probabalistic","Quick","Relativistic","Righteous","Runaway",
"Satellite","Secret","Sentient","System","Temporal","Vertical","Virtual","Water-cooled","Legendary","Elemental","Wireless","Argumentative","Mylenated","Agitated","Undercover",
"Spiritual","Depleted","Eigen","Marxist","Incendiary","Unnecessary","Sacrificial","Enlightened","Spatial"};
const char *nouns[46]=
{"Bitcoin","Butterfly Wings","Calculations","Capacitors","Comics","Computations","Cookies","Cpu","Estimate","Files","Gnomes","Hamsters","Harmonics","Heat Sink","Inversion Tables",
"Matrices","Metatables","Models","Motherboard","Ozone Layer","Paths","Pathways","Patriarchy","Phasers","Pizza","Polarity","Predictions","Procedures","Ram","Robot Uprising",
"Russians","Satellite Terrain Data","Skynet","Spirit Realm","Splines","Subsystems","System Processes","Time Dilation","Transform","Universe","Wavelets","Marsupials","Doorman",
"Power Cells","Human Suffering","Priors"};
const char *adverbs[27]=
{"Aggressively","Barely","Begrudgingly","Gently","Happily","Hastily","Quickly","Unenthusiastically","Timidly","Confidently","Imaginatively","Patiently","Thoroughly","Proficiently",
"Significantly","Roughly","Deliberately","Over-confidently","Majestically","Vivaciously","Vainly","Vaguely","Vacantly","Judgmentally","Frantically","Awkwardly","Carelessly"};




// ADDED: 06/06/2018
// adjust axis for slices acquired along an axis other than z
// TESTED: 6/7/2018 - Tested on X and Y acquired Simulated Data.
// Test Path: /share/dbp2123/dparker/Code/TestSTC
// Seems to work fine.
// X-TODO-X: Need to make sure that when the object it flipped, Slice order row 1 still corresponds to slice 1
// 06/08/18: Actually, I'm calling this done.  thinking about it, I did test it on data I flipped onto the x axis, and if this function
// didn't handle it correctly, then the timing file wouldn't have been correct, and the signal wouldn't have
// been perfectly lined up, as it was in the sim.  

void adjust_axis(volume4D<float>& timeseries, std::string axis, std::string stage){
  
  // We assume that the fmri is acquired along the z axis, with slices in the XY plane
  // if this isn't the case, we have to realign the matrix so that it is.
  // This can be memory/computationally intensive, and a better way might be just to make
  // three different correction loops, but that's a lot of typing and idk I don't really want to do it.
  // If there's a better way, please let me know.
  
  
  // If we're doing a forward transform, we're taking the data in and transforming it
  // from its original shape so the program can work on it.  If we're doing the reverse,
  // transforming it from our modified state back to its original
  // 6/25/18 - Tested on sequential and Multiband X and Y axis acquisition - PASSED
  
  
	if (stage=="forward")
	{
	  cout << "Adjusting for slice acquisition along "<< axis <<" axis"<< endl;
	  if (axis=="x"){
		timeseries.swapdimensions(3,2,-1);
	  }
	  else if (axis=="y"){
		timeseries.swapdimensions(1,3,2);
	  }
	  else if (axis=="z"){
		cout << "The default axis is already z you dummy.  Maximum effort."<< axis << endl;
	  }
	  else{
		cout << "Invalid Axis option for --axis:"<< axis << endl;
	  }
	}
	else if (stage=="reverse")
	{
	  cout << "Undoing Adjustment for slice acquisition along "<< axis <<" axis"<< endl;
	  if (axis=="x"){
		timeseries.swapdimensions(-3,2,1);
	  }
	  else if (axis=="y"){
		timeseries.swapdimensions(1,3,2);
	  }
	  else if (axis=="z"){
		cout << "The default axis is already z you dummy.  Maximum effort."<< axis << endl;
	  }
	  else{
		cout << "Invalid Axis option for --axis:"<< axis << endl;
	  }
	}

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
	  int ModSample=std::abs(shift);
	  
	  // If the shift if positive (shifting the signal to the right), then we want to DELAY the filter, add zeros to the END (Right hand side)		
	  
  
	  padd.assign(std::abs(shift),pFIR.back());

	  pFIR.insert(pFIR.end(),padd.begin(),padd.end());
	  
	  
	  if (shift<0)
	  {
		  std::reverse(pFIR.begin(),pFIR.end());
		  ModSample=0;
	  }
	   
	  
	  ColumnVector FIR_down_shift;
	  ColumnVector FIR_down;
	  FIR_down_shift.ReSize(firLen);
	  FIR_down.ReSize(firLen);
	  lenF=FIR_down_shift.Nrows();	
	  firLen=1;
	  
	  for (unsigned i = 0; i< SamplePoints.size(); i++)
	  {
		  FIR_down(firLen)=FIR->operator[](SamplePoints[i]);
		  SamplePoints[i]+=ModSample;
		  FIR_down_shift(firLen)=pFIR[SamplePoints[i]];
		  firLen+=1;
	  }
	  
	  // If the shift if negative (shifting the signal to the left), then we want to add the zeros to the beginning (flip the signal)		
  
	  
	  
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

  int mod(int a,int b) {
	int c = a % b;
	return (c < 0) ? c + b : c;
  }

void make_timeseries20hz(volume4D<float> *timeseries, volume4D<float> *timeseries20hz, std::vector<float> *FIR, float TR, float Hf, int padlen)
{
  //filter_timeseries20hz(&cattimeseries, &FIR, &timings20hz)
	
// 		9/27/16 - modified Kaiser Window resampling algorithm and convolution filtering routine.
// 		Now matches output from old code almost perfectly (10e-3 error)
		
// 		6/25/18 - spent the past week adding this function to output the filtered function in 20hz
//		6/25/18 - STC and 20hz output on sequential z-acquired data - PASSED
	
	
	int linepoint=0;
	int xx=timeseries->xsize();
	int yy=timeseries->ysize();
	int zz=timeseries->zsize();
	int lenT = timeseries->tsize();
	int lenF = FIR->size();
	int lenT20 = timeseries20hz->tsize();
	float skip = Hf*TR;
	int intskip=int(round(skip));
	float tspan = lenF/Hf;
	float firStart=1*tspan/2.0;
	int nlow=0;
	int nhigh=0;
	int ntotal=(int) floor(tspan/TR-1);
	float toffset=padlen*TR;
	int shift=0;
	  
	volume<float> currentsum(xx,yy,zz);
  
	  float t20=0.0;
	  int point=0;
	  int nt=0;

	  
	  for (int t20n=0;t20n<lenT20;t20n++)
	  {
		t20=t20n/Hf;
		shift=mod((int) t20n+(int)lenF/2.0,(int)(Hf*TR));
		
		nlow=(int) ceil((t20+toffset-firStart)/TR);
		nhigh=(int) floor((t20+toffset+firStart)/TR);
		ntotal=nhigh-nlow;
		
		currentsum=0.0;
		

		point=((int) round(lenF-shift))-1;
		nt=0;
		float FIRsum=0.0;
		

		while (point>=0)
		{
		  linepoint=nhigh-nt;

		  currentsum+=((FIR->operator[](point))*(timeseries->operator[](linepoint)));
		  nt++;
		  FIRsum=FIRsum+FIR->operator[](point);
		  point=point-intskip;
		  
		}
		timeseries20hz->operator[](t20n)=currentsum;
	  }
}


void output_message(){
	if (verbose.set())
	{
	  int size;
	  int v1;
	  std::string message="";
	  Matrix rmat=unifrnd(1,3,0,100);
	  Matrix choose;
	  int bins[10]={0,0,0,0,0,0,0,0,0,0};
	  int prob=0;
	  
	  prob=(int) floor(rmat(1,1));
	  if (prob<=20){
		 
		  size = *(&oneoff + 1) - oneoff;
		  choose=unifrnd(1,1,0,size-1);
		  v1=choose(1,1);
		  message=message+oneoff[v1]+" ";
		  //cout<<message<<endl;
		  
	  }
	  
	  else{
		  prob=(int) floor(rmat(1,2));
		  if (prob<=10){
			  size = *(&adverbs + 1) - adverbs;
			  choose=unifrnd(1,1,0,size-1);
			  v1=choose(1,1);
			  message=message+adverbs[v1]+" "; 
		  }
	
		  size = *(&verbs + 1) - verbs;
		  choose=unifrnd(1,1,0,size-1);
		  v1=choose(1,1);
		  message=message+verbs[v1]+" ";
	
		  
		  size = *(&adjectives + 1) - adjectives;
		  choose=unifrnd(1,2,0,size-1);
		  v1=choose(1,1);
		  message=message+adjectives[v1]+" ";
		  
		  prob=(int) floor(rmat(1,3));
		  if (prob<=20){
			  v1=choose(1,2);
			  message=message+adjectives[v1]+" ";
		  }        
		  
		  size = *(&nouns + 1) - nouns;
		  choose=unifrnd(1,1,0,size-1);
		  v1=choose(1,1);
		  message=message+nouns[v1]+" "; 
  
	  }
	  cout<<message<<endl;
  }
}



void make_timings(Matrix *timings, Matrix *orders, int zs)
{
	
	// This is  the timing file case
	if ( timing.set() )
	{
		
		if ( refslice.set()||reftime.set() )
		{
			std::cout<<"When using a Slice Timing file, the times in the file supercede all other settings.  The reference slice/time specified will be ignored"<<std::endl;
			std::cout<<"If you wish to align data to a specific slice, please make that adjustment in the Slice Timing File, or omit the slice timing file."<<std::endl;
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
	
	// this is the order file case
	else if ( order.set() )
	{
		//Slice Order File, Default Reference Slice Tested 9/23/16 - Passed
		//Slice Order File, Custom Reference Slice Tested 9/23 - Passed		
		
		// 9/23/16 - Will not Run With Multiband, only slice timing file (User Burden, Deal With It)
		
		//Slice Order File, Default Reference Slice Tested 10/10/16 - Passed
		//Slice Order File, Custom Reference Slice Tested 10/10/16 - Passed
		
		// 5/12/17 - Working on glitch in sequential and even odd interleaev (shift times not working)
		// - completed, Added "else" case that catches when no ref slice or ref time is provided
		
		// 6/6/2018 - working on allowing multiband to be used with slice order.
		// OK here's the problem.  the slice order file is NOT a file that tells you what order a slice was acquired in.
		// It tells you what the order in which the slices were acquired.  Confused?  Let's take a look:
		// Note: slice order files are read from top to bottom
		//
		// This is the typical slice order format.  It means slice "val" was acquired "row"
		//
		//   Row:  Val:
		//   1     1
		//   2     3
		//   3     5
		//   4     7  
		//   5     2
		//   6     4
		//   7     6
		//   8     8
		//
		//  In this example, "row" indicates the order in which the slices are acquired, "val" indicates the slice number
		// SO, row 2 val 3 means slice 3 was acquired 2nd.  row 5 val 2 means slice 2 was acquired 5th, and so on
		//
		//  Multiband would need the following format, which means slice "row" was acquired "val"
		//
		//   Row:  Val:
		//   1     1
		//   2     3
		//   3     2
		//   4     4  
		//   5     1
		//   6     3
		//   7     2
		//   8     4
		//
		// Now in this case, row 2 val 3 means slice 2 was acquired 3rd.  Row 5 val 1 means slice 5 was acquired 1st
		// So now I have to make a smart thing to detect if it's multiband. THANKS A LOT YOU JERKS
		
		
		int tmx;
		float shift;
		bool multiband=false;
		
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

		if (tmx<zs){
			std::cout<<"Multiband Mode Detected"<<std::endl;
			output_message();
			multiband=true;
		}
		
		
		shift=TR.value()*1.0/(tmx+1);
		
		
		// create a time list - the time of the ith acquisition, one value for all z's		
		Matrix TimeList;
		TimeList.ReSize(zs,1);
		
		for ( int i=1; i<=zs; i++ )
		{
			TimeList(i,1)=(float) (i-1)*shift;
		}

		// 9/23/16 - Made This loop more efficient
		// so normally, slice order is: row 2 means this slice was acquired 2nd, row 3 means third...etc
		// BUT the timings file is row 1 is the time of slice 1 acquisiton, row 2 is slice 2 time, etc
		// so order[3]= slice acquired 3rf, timings[order[3]] is how we index that slice, and
		// timelist[3] is the time of whichever slice was acquired 3rd.
		
		if (multiband)
		{
		  // if we've detected multiband, then we assume that slice order is now "row 2 val 3 means slice 2 was acquired 3rd"
		  // so order[3] means slice 3 is acquired VAL
		  // timing [3] is timing of slice 3
		  // TimeList[order[3]] = time of slice 3
		  for (int j=1;j<=zs;j++)
		  {
			timings->operator()(j,1)=TimeList(orders->operator()(j,1),1);
		  }
		}
		else
		{
		  // otherwuse use original algorithm
		  
		  for (int j=1;j<=zs;j++)
		  {			
			  timings->operator()(orders->operator()(j,1),1)=TimeList(j,1);					
		  }
		}
		
		// this is reference slice stuff, it shoudl work for multiband, but:
		// TODO: Make sure ref doesn't fall outside of TR
		if ( refslice.set() )
		{		
			timings->operator-=(timings->operator()(refslice.value(),1));
		}
		else if ( reftime.set() )
		{
			timings->operator-=(reftime.value());
		}
		
	}
	
	// this is the case if there's just a reference slice given (TR and itl assumed provided)
	else if (refslice.set() )
	{
		// Create Timing File Tested 9/22/16 - Succesful
		// With Reference Tested 9/22/16 - Succesful
		
		float dt;
		Matrix IntSeq;
		Matrix TimeList;
		IntSeq.ReSize(zs,1);		
		TimeList.ReSize(zs,1);
		dt=TR.value()/(zs+1);
		int counter;
		counter=1;

		output_message();
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
		TimeList-=TimeList(refslice.value(),1);
		*timings=TimeList;
		

	}
	
	// this is the case if there's just a reference time given (TR and itl assumed provided)
	else if (reftime.set() )
	{
		// 11/01/16 - added, need to test.
		
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
		TimeList-=reftime.value();
		*timings=TimeList;
		
    // 5/12/17 - added to catch when no ref slice or ref time is provided (assumes slice 1 for reference)
	}
	else
	{
		float dt;
		Matrix IntSeq;
		Matrix TimeList;
		IntSeq.ReSize(zs,1);		
		TimeList.ReSize(zs,1);
		dt=TR.value()/(zs+1);
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
		TimeList-=TimeList(1,1);
		*timings=TimeList; 

	  
	  
	}
	
	
}

std::vector<int> Padd_Timeseries(ColumnVector *voxeltimeseries, ColumnVector *cattimeseries, int no_volumes)
{
  //std::cout<<"flippingTS"<<std::endl;
	ColumnVector fliptimeseries = voxeltimeseries->Reverse();
	fliptimeseries=fliptimeseries.Rows(2,no_volumes-1);
	std::vector<int> cutoff(2);
	int cutLeft;
	int cutRight;
	
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
	
	// 11/28/16  Maybe this "+2" has something to do with the negative bold problem...?
	// No, just checked.  
	
	cutLeft = lents-rangelh+2;
	cutRight = cutLeft+no_volumes-1;
	
	if (cutRight-cutLeft != no_volumes-1)
	{
		cout<<"Problem in padding, check output carefully"<<endl;
	}
	
	
	//std::cout<<"ConcateTimeseries"<<std::endl;
	*cattimeseries=fliptimeseries.Rows(rangelh,lents)&*voxeltimeseries&fliptimeseries.Rows(1,rangerh);	
	cutoff[0]=cutLeft;
	cutoff[1]=cutRight;
	return cutoff;
}

int shift_volume()
{
	
	// 10-10-16 Added check for some combination of optional arguments
	// Now you idiots can't mess up your data by leaving out these settings.
	// ...but I'm sure you'll still manage to find a way >:(
	
	if (!TR.set()&&!order.set()&&!timing.set()&&!lpf.set()&&!hpf.set())
	{
		std::cout<<"Must Specify either a TR, a slice Order file,a slice timing file,\n or indicate that you only wish to preform filtering"<<std::endl;
		return -1;
	}
	
	
	Matrix timings;
	Matrix orders;
	std::vector<float> FIR;
	

	int no_volumes = 0; 
	int xx = 0;
	int yy = 0;
	int zz = 0;
	int HRhi=20;
	volume4D<float> timeseries;
	

	  //volume4D<float> timeseries;
	    if (input.set())
		{
		  
		  if (true) { cout << "Reading input volume" << endl; }  // DO NOT MESS WITH THIS IF STATEMENT ITS VERY IMPORTANT
		  output_message();
		  read_volume4D(timeseries,input.value());
		  
		  
		  
		  // If we set which axis we're using, correct for it so the triple for loops will work
		  if ( axis.set() ){
			// ADDED 06/06/2018
			// Tested
			adjust_axis(timeseries,axis.value(),"forward");
		  }
		no_volumes = timeseries.tsize(); 
		xx = timeseries.xsize();
		yy = timeseries.ysize();
		zz = timeseries.zsize();
		  
		} else if (out.set()) {
		  cerr << "Must specify an input volume (-i or --in) to generate corrected data." << endl;
		  return -1;
		}
	  



	std::cout<<"Create timing Arrays"<<std::endl;
	output_message();
	timings.ReSize(timeseries.zsize(),1);
	orders.ReSize(timeseries.zsize(),1);	
	make_timings(&timings,&orders,timeseries.zsize());
	
	string Output=out.value();
	string directory=Output.substr(0,Output.find_last_of('/')+1);

	float cutoff=cf.value();
	
	float samplingrate=(float) (zz/TR.value());
	

	
	float stopgain=-28;
	double transwidth=.1;
	int PassZero=1;
	int skip=samplingrate*TR.value();
	

	// 01/16/17 - HPF can't operate the same way LPF does (on zero-padded data), so just filter normally.
	// 01/16/17 - Testing HPF operation now...
	if (hpf.set())
	{
		PassZero=0;
		samplingrate=(float) (1.0/TR.value());
		skip=1;
	}
	
	//std::cout<<"Pass Zero: "<<PassZero<<std::endl;
	std::cout<<"Generate Filter START"<<std::endl;
	output_message();
	window kaiser(cutoff,samplingrate,stopgain,transwidth,PassZero,TR.value());
	FIR=kaiser.get_fir();
	std::cout<<"Generate Filter FINISHED - success\n"<<std::endl;
	output_message();
	
	// I think this is just initializing the values that will be used in the loop
	ColumnVector voxeltimeseries = timeseries.voxelts(1,1,1);
	ColumnVector cattimeseries;
	std::vector<int> padcut(2);
	int cutLeft;
	int cutRight;
	
	float span;
	float mn;
	float mn2;
	float span2;
	
	// 10/19/16 - "lpf" option now specifies JUST filtering.  This part sets
	// the slice delay to be zero for all of them (No shifting.  NO SHIFTING!)
	
	if (lpf.set()||hpf.set())
	{
		for (int tm=1;tm<=zz;tm++)
		{
			timings(tm,1)=0;
		}
	}
	write_ascii_matrix(directory+"TimingFile.txt", timings, 6);
	std::cout<<"Timing file wrote to: "<<directory<<"TimingFile.txt\n"<<std::endl;
	std::cout<<"Filtering START..."<<std::endl;
	output_message();
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
				
				// 01/16/17 - re-added original mean so HPF functions correctly...but there seems to be an offset in the final data.
				voxeltimeseries+=(mn);
				
				// 9/27/16 - Filter changes mean and span - added code to maintain mean and span.
				// Also checks to see if the span is zero.  If so, do not filter.
				
				if ( span!=0 )
				{
				  
					if (hpf.set())
					{
					  //std::cout<<"StartingPadd"<<std::endl;
					  padcut=Padd_Timeseries(&voxeltimeseries, &cattimeseries,no_volumes);
					  //std::cout<<"Done"<<std::endl;
					  cutLeft=padcut[0];
					  cutRight=padcut[1];
					  
					  //std::cout<<"Writing File"<<std::endl;				
					  //write_ascii_matrix(cattimeseries,"/home/dparker/Desktop/MyOutput/FiltershiftTest/testTS.txt", 16);
					  //std::cout<<"WroteFile"<<std::endl;
					  //std::exit(0);		
					  filter_timeseries(&cattimeseries, &FIR, (int)floor(timings(slice,1)*(float)samplingrate+0.5),skip,slice);
					  cattimeseries=cattimeseries.Rows(cutLeft,cutRight);					  
					  
					  
					  
					  
					  
					  
					}
					//else
					//{
					//  fliptimeseries = voxeltimeseries.Reverse();
					//  fliptimeseries=fliptimeseries.Rows(2,no_volumes-1);				
					//  cattimeseries=fliptimeseries.Rows(rangelh,lents)&voxeltimeseries&fliptimeseries.Rows(1,rangerh);
					//  
					//  //std::cout<<"Writing File"<<std::endl;				
					//  //write_ascii_matrix(cattimeseries,"/home/dparker/Desktop/MyOutput/FiltershiftTest/testTS.txt", 16);
					//  //std::cout<<"WroteFile"<<std::endl;
					//  //std::exit(0);		
					//  filter_timeseries(&cattimeseries, &FIR, (int)floor(timings(slice,1)*(float)samplingrate+0.5),skip,slice);
					//  cattimeseries=cattimeseries.Rows(cutLeft,cutRight);
					//}
					//
  
					// 01/16/17 - Only remean and adjust span for LPF.
					if (!hpf.set())
					{
					  mn2=cattimeseries.Sum()/cattimeseries.Nrows();
					  cattimeseries-=(mn2);
					  span2=cattimeseries.Maximum()-cattimeseries.Minimum();
					  cattimeseries/=span2;
					  cattimeseries*=span;
					  cattimeseries+=mn;
					}					
					//std::cout<<"resetting TS"<<std::endl;
					//std::cout<<cattimeseries.Nrows()<<std::endl;
					timeseries.setvoxelts(cattimeseries,x_pos,y_pos,slice-1);
					//std::cout<<"Done"<<std::endl;
				}
			}
		}
	}
	
	
	  // If we set which axis we're using, un-correct for it sso the volume looks the same as the input
	if ( axis.set() ){
	  // ADDED 06/06/2018
	  // Tested
	  adjust_axis(timeseries,axis.value(),"reverse");
	}
	
	
	
	std::cout<<"Filtering FINISHED - success\n"<<std::endl;
	output_message();
	std::cout<<"Writing Output Volume"<<std::endl;
	output_message();
	write_volume4D(timeseries,out.value());	
	

//######################################################################################################
// input stuff incase Ray wants 20Hz.  I'm not clever enough to do this a different way.
// 06/13/2018 - Adding hires part.  
//######################################################################################################
	if (hires.set())
	{
	  
	  std::vector<float> FIR20;
	  
		int orig_t = timeseries.tsize();
		xx = timeseries.xsize();
		yy = timeseries.ysize();
		zz = timeseries.zsize();
		
		volume4D<float> timeseries20hz(xx,yy,zz,round((orig_t+1)*TR.value()*HRhi));		

		no_volumes = timeseries20hz.tsize();

		int skip=(int) round(TR.value()*HRhi);
		float timingshift=1.0/HRhi;
		
		cutoff=1.0/(2*TR.value());
		
		
		window kaiser20hz(cutoff,HRhi,stopgain,transwidth,PassZero,TR.value());
		FIR20=kaiser20hz.get_fir();
		float startval=FIR20[0];
		//kaiser20hz.print_info();
		std::cout<<"Generate Filter FINISHED - success\n"<<std::endl;
		output_message();
		int lenFIR=FIR20.size();
		int padlen=ceil(lenFIR/(HRhi*TR.value()));
		std::cout<<"Padlen:\t"<<padlen<<std::endl;
		volume4D<float> padded_timeseries(xx,yy,zz,orig_t+padlen*2);
		
		for (int i=0; i<padlen; i++)
		{
		
		  padded_timeseries[i]=timeseries[padlen-i];
		  padded_timeseries[orig_t+padlen+i]=timeseries[orig_t-i-2];
		}
		
		for (int i=padlen;i<padlen+orig_t;i++)
		{
		
		  padded_timeseries[i]=timeseries[i-padlen];
		}
		for (int i=0;i<lenFIR;i++){
		  FIR20[i]=FIR20[i]-startval;
		}
		padded_timeseries.swapdimensions(-1,2,3);  // for some reason the x axis got flipped here when I assigned the TS to padded_timeseries
		// possibly a problem with FSL's assignment code, possibly a problem with my hacked to gether code.  
		//std::for_each(FIR20.begin(),FIR20.end(),[](float& d ) {d-=startval;});

	  make_timeseries20hz(&padded_timeseries, &timeseries20hz, &FIR20, TR.value(),HRhi, padlen);

	  std::cout<<"Writing Highres Output Volume"<<std::endl;
	  output_message();
	
	  write_volume4D(timeseries20hz,out.value()+"_Highres");	
	  
	}

  return 0;
}





int main (int argc,char** argv)
{
  
  srand(seed);
  
  output_message();
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
	options.add(refslice);
	options.add(reftime);
	options.add(lpf);
	options.add(hpf);
	options.add(axis);
	options.add(hires);
	options.add(verbose);
	
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
  
  // 11/01/16 - Added catch to prefent cutoff being larger than nyquist.  For some reason,
  // This would really mess up the filtering...
  if (cf.value()==0||cf.value()>(1.0/(2.0*TR.value())))
  {
	std::cout<<"Adjusting Cutoff to be <= Nyquist"<<std::endl;
	float cutoff=(float) 1.0/(2.0*TR.value());
	std::ostringstream ss;
	ss << cutoff;
	std::string s(ss.str());
	cf.set_value(s);
  }
  
  if (reftime.value()>TR.value())
  {
	cerr<<"Reference time set by --rt must fall within 0 <= rt < TR"<<std::endl;
	exit(EXIT_FAILURE);
	
  }
  output_message();
  std::cout<<"Processing START"<<std::endl;
  int retval = shift_volume();
  std::cout<<"Processing FINISH - success\n"<<std::endl;
  
  if (retval!=0) {
	cerr << endl << endl << "Error detected: try -h for help" << endl;
  }
	
	
  return retval;
}
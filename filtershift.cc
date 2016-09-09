#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>

#include "miscmaths/optimise.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "utils/options.h"
#include "miscmaths/kernel.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace Utilities;

string title="FilterShift \nQNL's tool for optimal slice timing correction";
string examples="filtershift --in <InputFile> --tr <TR> [options]";

Option<string> input(string("--in"), string(""),
			  string("\tfilename the input image to perform STC on\n"),
			  true, requires_argument);

Option<float> TR(string("--TR"), 2.0,
			  string("\tSet the TR of the original fMRI data in seconds\n"),
			  true, requires_argument);

Option<int> interleave(string("-i,--interleave"), 0, 
		     string("\tset the interleave parameter, or how many slices are incramented between acquisitions\
                    \n\t\t\t\t\t1 = sequential acquisition\
                    \n\t\t\t\t\t2 = even/odd acquisition (acquire every second slice)\
                    \n\t\t\t\t\t3 = acquire every third slice, etc\n"), 
		     false, requires_argument);

Option<string> out(string("-o,--out"), string(""), 
		     string("\tSpecify an output file name - all working directories will be created in the parent directory specified here.\
                    \n\t\t\t\t\tLeave blank to run in the parent directory of <InputFile>\n"), 
		     false, requires_argument);

Option<int> start(string("-s,--start"), 0,
			  string("\tSet the starting slice - The slice that was acquired first in the sequence.  Default is slice 1, the bottom most slice.\
                     \n\t\t\t\t\tThis starts the interleave from that slice.  If your interleave parameter is '1'\
                     \n\t\t\t\t\tand your starting slice is '3', your slice acquisition sequence will be modled as:\
                     \n\t\t\t\t\t3\n\t\t\t\t\t5\n\t\t\t\t\t7\n\t\t\t\t\t9...\n"),
			  false, requires_argument);
              
Option<int> direction(string("-d,--direction"), 1, 
		     string("\tvalue 1 or -1.  Set the direction of slice acquisition.\
                    \n\t\t\t\t\t 1: implies ascending slice acquisition from starting slice: (3,5,7,9...)\
                    \n\t\t\t\t\t-1: implies descending slice acquisition from starting slice: (7,5,3,...)\n"), 
		     false, requires_argument);

Option<string> order(string("--order"), string(""),
			 string("\t\tSlice Order File.  This file is the order in which each slice was acquired. each row represents the order in which that slice was acquired.\
                    \n\t\t\t\t\t For example, '1' in the first row means that slice 1 was aquired first. '20' in the second row means that slice 2 was acquired 20th.\
                    \n\t\t\t\t\t If present, all interelave parameters are ignored, and slices are shifted using the slice order file\
                    \n\t\t\t\t\t we refer to the bottom slice in the image as slice 1, not slice 0\n"),
			 false, requires_argument);

Option<string> timing(string("--timing"), string(""),
			 string("\tSlice Timing File.  This file is the time at which each slice was acquired relative to the first slice. each row represents the time at which that slice was acquired.\
\n\t\t\t\t\t For example, '0' in the first row means that slice 1 was aquired first, and will be shifted 0 seconds.  '0.5' in the second row\
                      \n\t\t\t\t\t means that slice 2 was acquired 0.5 seconds after the first slice, and will be shifted 0.5s.\
                      \n\t\t\t\t\t If present, all interelave parameters are ignored, and slices are shifted using the slice timing file.\
                      \n\t\t\t\t\t If you put both a slice timing and a slice order, the program will yell at you and refuse to run.\n"),
			 false, requires_argument);

Option<int> ref(string("-r,--reference"), 0,
			 string("\tSet the Reference slice\n\tThis is the slice the data is aligned to.  Default is the first slice\n"),
			 false, requires_argument);
              
Option<int> cores(string("-c,--cores"), 1,
			 string("\tthe number of cores free on your computer for parallel processing.  The default is one less than the number of cores your computer has.\
                    \n\t\t\tIt is highly recommended that you use parallel processing to speed up this operation\n"),
			 false, requires_argument);

Option<float> mem(string("-m,--mem"), 4.0,
			 string("\tamount of free memory you have, in GB, so the system knows how much load it can give each core. The default is 4\n"),
			 false, requires_argument);
              
Option<bool> force(string("--force"), false,
			 string("\tForces the program to proceed even if it estimates a memory error, or if it detects that STC has already been run.  Use -force if you need to re-run with new parameters\
                    \n\t\t\t\t\tWARNING: THIS CAN POTENTIALLY CRASH YOUR COMPUTER\n\n\n"),
			 false, no_argument);              
            
            
  Option<bool> help(string("-h,--help"), false,
		  string("\tdisplay this message\n"),
		  false, no_argument);

int shift_volume()
{
	cout << "Testing A Load and swap";
  volume4D<float> timeseries;
  

  if (input.set()) {
    if (true) { cout << "Reading input volume" << endl; }
    read_volume4D(timeseries,input.value());
    if (!out.set())
      out.set_value(input.value() + "_st");
  } else if (out.set()) {
    cerr << "Must specify an input volume (-i or --in) to generate corrected data." 
	 << endl;
    return -1;
  }

  int no_volumes = timeseries.tsize();
  print_info(timeseries,input.value());
  Matrix imgmat;
  imgmat=timeseries.matrix();
  int NR,NC;
  NR=imgmat.Nrows();
  NC=imgmat.Ncols();
  cout << NR << endl;
  cout << NC << endl;
  cout << imgmat(50,100) << endl;
  return 0;
}





int main (int argc,char** argv)
{
  //Tracer tr("main");
  
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
    options.add(timing);
    options.add(ref);
    options.add(cores);
    options.add(mem);
    options.add(force);

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

  int retval = shift_volume();
  
  if (retval!=0) {
    cerr << endl << endl << "Error detected: try -h for help" << endl;
  }
	
	
  return retval;
}
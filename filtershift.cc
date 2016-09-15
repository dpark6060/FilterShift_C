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
#include "Window.h"

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

Option<int> ref(string("-r,--reference"), 1,
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





void make_timings(Matrix *timings, Matrix *orders, int zs)
{
	if ( timing.set() )
	{

		Matrix TempTimings;
		TempTimings.ReSize(zs,1);
		float tmn;
		float tmx;		
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
		
		std::cout<<TempTimings<<std::endl;
		tmx=TempTimings.Maximum();
		
		
		for ( int i=1; i<=zs; i++ )
		{
			tmn=Minimum(TempTimings);

			if ( tmn >tmx )
			{
				break;
			}
			
			for ( int j=1;j<=zs; j++ )
			{
				std::cout<<"timings("<<j<<") = "<<timings->operator()(j,1)<<std::endl;
				
				if ((float) timings->operator()(j,1)==tmn)
				{
					orders->operator()(j,1)=i;					
					TempTimings(j,1)=tmx*2;
					std::cout<<TempTimings<<std::endl;
				}
			}
		}
		
		
		
		
	}	
	else if ( order.set() )
	{
		std::cout<<"Order File"<<std::endl;
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
		
		
		orders->operator-=(1);
		tmx=orders->Maximum();
		shift=(TR.value()*1.0/zs);
		for ( int i=1; i<=zs; i++ )
		{
			timings->operator()(i,1)=orders->operator()(i,1)*shift;
		}
		
		std::cout<<ref.value()<<std::endl;
		timings->operator-=(timings->operator()(ref.value(),1));
		
	}
	else
	{
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
		
		std::cout<<"Made"<<std::endl;

	}
	
	
	
}


int shift_volume()
{
	volume4D<float> timeseries;
	Matrix timings;
	Matrix orders;
	
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
  
	std::cout<<"Read Succesful"<<std::endl;
	timings.ReSize(timeseries.zsize(),1);
	orders.ReSize(timeseries.zsize(),1);	
	make_timings(&timings,&orders,timeseries.zsize());
	
	//std::cout<<"TimingFile:"<<std::endl;
	//for (int i=1;i<=timeseries.zsize();i++)
	//{
	//	std::cout<<timings(i,1)<<std::endl;
	//}
	//std::cout<<"OrderFile:"<<std::endl;
	//for (int i=1;i<=timeseries.zsize();i++)
	//{
	//	std::cout<<orders(i,1)<<std::endl;
	//}
	
	int no_volumes = timeseries.tsize();
	int xx = timeseries.xsize();
	int yy = timeseries.ysize();
	int zz = timeseries.zsize();
	std::cout<<xx<<std::endl;
	volume4D<float> filterseries(xx,yy,zz,no_volumes*2-2);
	
	
	
	// Test Window Function
	float cutoff=0.2;
	int samplingrate=20;
	float stopgain=-60;
	float transwidth=.08;
	
	//window::window kaiser(cutoff,samplingrate,stopgain,transwidth);
	//kaiser.print_info();
	//kaiser.get_fir();
	
	
	for (int slice=1; slice<=zz; slice++) {
	
		for (int x_pos = 0; x_pos < xx; x_pos++)
		{
			for (int y_pos = 0; y_pos < yy; y_pos++)
			{
		
				ColumnVector voxeltimeseries = timeseries.voxelts(x_pos,y_pos,slice-1);
				ColumnVector fliptimeseries=voxeltimeseries.Reverse();
				
				//std::cout<<"\n\n"<<std::endl;
				//std::cout<<voxeltimeseries(199)<<std::endl;
				//std::cout<<voxeltimeseries(200)<<std::endl;
				
				fliptimeseries=fliptimeseries.Rows(2,no_volumes-1);
				
				int lents=fliptimeseries.Nrows();
				int rangelh;
				int rangerh;
				
				if ( lents % 2 == 0 )
				{
					rangelh=no_volumes-lents/2-2;
					rangerh=lents/2-1;
					
					
				}
				else
				{
					rangelh=no_volumes-floor(lents/2)-2;
					rangerh=ceil(lents/2);
				}
				
				//std::cout<<"Catting File"<<std::endl;
				//
				//std::cout << "rangelh: " << rangelh << std::endl;
				//std::cout << "no_volumes: " << no_volumes << std::endl;
				//std::cout << "rangerh: " << rangerh << std::endl;
				//std::cout << "fliptimeseres.Nrows: " << fliptimeseries.Nrows() << std::endl;
				//std::cout << "rangelh: " << rangelh << std::endl;
				//
				//
				//ColumnVector cattimeseries=fliptimeseries.Rows(rangelh,no_volumes-2)&voxeltimeseries&fliptimeseries.Rows(1,rangerh);
				//
				//std::cout << "Len cattimeseries: " << cattimeseries.Nrows() << std::endl;
				//std::cout << "Len filterseries: " << filterseries.tsize() << std::endl;
				//std::cout << "mint filterseries: " << filterseries.mint() << std::endl;
				//std::cout << "maxt filterseries: " << filterseries.maxt() << std::endl;
				
				//std::cout<<"Writing File"<<std::endl;
				//
				//write_ascii_matrix(cattimeseries,"/home/dparker/Desktop/MyOutput/FiltershiftTest/testTS.txt", 16);
				//std::cout<<"WroteFile"<<std::endl;
				//return 1;
				filterseries.setvoxelts(cattimeseries,x_pos,y_pos,slice-1);
				
				
				
				//std::cout<<voxeltimeseries.Nrows()<<std::endl;
				//std::cout<<voxeltimeseries(197)<<std::endl;
				//std::cout<<voxeltimeseries(198)<<std::endl;
				

				
				//for (int time_step=1; time_step <= no_volumes; time_step++)
				//{
				//	// interpseries(time_step) = interpolate_1d(voxeltimeseries, time_step - offset);
				//	//interpseries(time_step) = kernelinterpolation_1d(voxeltimeseries, time_step + offset, userkernel, 7);
				//	//std::cout<<"here"<<std::endl;
				//}
				//timeseries.setvoxelts(interpseries,x_pos,y_pos,slice-1);
			}
		
			//if (verbose.value())
			//cerr << "Slice " << slice << " offset " << offset << endl;
		}
	}
	
	write_volume4D(filterseries,"/home/dparker/Desktop/MyOutput/FiltershiftTest/FiltSeries.nii");
		//std::cout<<flipseries(1,1)<<std::endl;
		//std::cout<<flipseries(1,1000)<<std::endl;
		
		
		//print_info(timeseries,input.value());
	
	
	
		//cout << imgmat(50,100) << endl;
		

  
  
  
  
  
  
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
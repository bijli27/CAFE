#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)
#include <stdio.h>
#include "caHdl.h"
#include <time.h>
#include <ctime>
#include <iostream>
#include <iomanip>

//###LB: Unless indicated all angles are calculated in degrees.

#define SOLCAPARVERSION 1.0		//version of the program

int main( int argc, char *argv[] )
{

	std::time_t time1 = std::time(nullptr);
	std::cout << "Start Time is "
		<< std::asctime(std::localtime(&time1));

#ifdef _MEMLEAK_
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif
	ostringstream message;// ostring class from standard library
	caHdlP ca = new caHdl(); // (not clear) caHdlP is a pointer of class caHdl

	/*if (argc < 2) 
		ca->err->reportError(ERRTXT("Not enough input arguments"));*/

	/*FILE *test;
	test = fopen("thisIsATestFile_parallel.txt","w");
	fclose(test);*/
	ca->initializeSolCA("InputFile_MultiLayer_750_withoutRotation.uds", 0 ); // calling the function  initializeSolCA with argument under class caHdl. object is ca and member is initializeSolCA of class caHdlp
	//ca->initializeSolCA("InputFile_SingleLayer_750.uds", 0);
	std::time_t time2 = std::time(nullptr);
	std::cout << "Start Parallel Process Time is "
		<< std::asctime(std::localtime(&time2));

	ca->startParallelProcesses(argc, argv);// (accessing member function in class ca dreferrecning of pointer ca nad then the grab the member of  class. it is as like (*ca).startParallelProcesses
	ca->prepareAutomata();// (accessing member function in  class ca)
	ca->startSimulation();// (accessing member function in  class ca)
	//ca->endParallelProcesses();
	//delete ca;

//	ca->endParallelProcesses();
#ifdef _MEMLEAK_
	_CrtDumpMemoryLeaks();
	_CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

	std::time_t time3 = std::time(nullptr);
	std::cout << "End Time is "
		<< std::asctime(std::localtime(&time3));

	std::cout << "Program Complete" << std::endl;
	return 0;
}
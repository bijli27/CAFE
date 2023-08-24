#include "applic.h"
#include <stdio.h>
#include "parallel.h"

void exitus(const char *s)   //exit function for abnormal termination
{
	printf("%s\n", s);
	MPI_Abort(MPI_COMM_WORLD, 1);  //calling abort method of MPI
	exit(1);
}

bool trigger::now(void)  		//trigger function for time dependent events
{
	if (never) return NOTYET;
	whenOut += 1;
	if (whenOut < frequency) 		return NOTYET;
	whenOut = 0;
	return NOW;
}

trigger::~trigger(void) 		//destructor of trigger class
{
}


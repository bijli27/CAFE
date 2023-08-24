#pragma once
#ifndef _simhdl_h_
#define _simhdl_h_

#include "io.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>

class simHdl : public io				// simHdl inherits io
{
public:
	simHdl( void );						//constructor of simHdl class
	~simHdl( void );					//destructor of simHdl class
	void endProgram(void);
};

#endif
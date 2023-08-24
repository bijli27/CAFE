#include "io.h"
#include "caHdl.h"


tempCellPool::tempCellPool(caHdl* own)
{
	this->owner = own;
	this->tcellraw = NULL;
	this->nodes = new tempNodePool(this);
	PNBucket.TCparallelNode = (tempCellP*)calloc(MIN_TCELLS_PERSUBDOMAIN, sizeof(tempCellP));
	PNBucket.size = 0;
	PNBucket.maxLen = MIN_TCELLS_PERSUBDOMAIN;
	timeEncodeFactor = 0;
	this->first=NULL;
	this->last=NULL;
	this->all=NULL;
}

tempCellPool::~tempCellPool(void)
{
	for(long i=0;i<this->ntcells;i++)
	{
		delete this->all[i];
		this->all[i] = NULL;
	}

	delete all;

	if(T_timer.timeSeries)
	{
		free(T_timer.timeSeries);
		T_timer.timeSeries = NULL;
	}
	if(xdiscr)
	{
		free(xdiscr);
		xdiscr = NULL;
	}
	if(ydiscr)
	{
		free(ydiscr);
		ydiscr = NULL;
	}
	if(zdiscr)
	{
		free(zdiscr);
		zdiscr = NULL;
	}
	if(PNBucket.TCparallelNode)
	{
		free(PNBucket.TCparallelNode);
		PNBucket.TCparallelNode=NULL;
	}
	if(nodes)
	{
		delete nodes;
		nodes = NULL;
	}
}

Real tempCellPool::getCellNucleationTime(cellP c)
{
	QUICKASSERT(((uint32_t)c) & 0x01);
	/*Real z = ((uint32_t)c) & 0x01;
	if (z ==0)
	{
		uint32_t nucTime = (((uint32_t)c) & 0x7FFE) >> 1;
	}*/

	uint32_t nucTime = (((uint32_t)c) & 0x7FFE) >> 1;
	return ((Real)nucTime) / timeEncodeFactor;
}

Real tempCellPool::getCellNucleationTime(uint_fast64_t data)
{
	uint_fast64_t dataFiltered = data & 0x1FFFE0;
	dataFiltered >>= 5;

	return ((Real)dataFiltered) / timeEncodeFactor;
}

void tempCellPool::setTimeEncodeFactor(void)
{
	long len = T_timer.size - 1;
	Real maxTime = T_timer.timeSeries[len];
	maxTime = maxTime * 2;

	timeEncodeFactor = (Real)TIMEDISCRETIZATION / maxTime; //TIMEDISCRETIZATION is 16383. we are dividing the number of time steps with the maximum time.Why??
}

void tempCellPool::addTempCelltoBucket(tempCellP Tcell)
{
	QUICKASSERT(this->PNBucket.TCparallelNode);
	
	int i = 0;

	for (i = 0; i < PNBucket.size; i++)
		if (PNBucket.TCparallelNode[i] == Tcell) return;

	if (i < PNBucket.maxLen)
	{
		PNBucket.TCparallelNode[i] = Tcell;
		PNBucket.size++;
		return;
	}

	size_t newsize = (size_t)PNBucket.maxLen + (size_t)MIN_TCELLS_PERSUBDOMAIN;
	tempCellP* newMemory = (tempCellP*)calloc(newsize, sizeof(tempCellP));
	tempCellP* oldMemory = PNBucket.TCparallelNode;

	memcpy(newMemory, oldMemory, (size_t)PNBucket.maxLen * sizeof(tempCellP));

	PNBucket.maxLen = (short)newsize;
	PNBucket.TCparallelNode = newMemory;
	
	PNBucket.TCparallelNode[i] = Tcell;
	PNBucket.size++;
	
	free(oldMemory);
}

long tempCellPool::addTempCell(tempCellP Tc)
{
	return 0;
}

void tempCellPool::setTtimer(Real minTime)
{
	long len = T_timer.size;

	for (long i = 0; i < len-1; i++)
	{
		if (minTime >= T_timer.timeSeries[i] && minTime <= T_timer.timeSeries[i + 1])
		{
			T_timer.key = i;
			return;
		}
	}
}

temperatureNode::~temperatureNode(void)
{
	if(Temperature)
	{
		free(Temperature);
		Temperature = NULL;
	}
}

temperatureNode::temperatureNode(Real x, Real y, Real z, int identity)
{
	pos.x = x;
	pos.y = y;
	pos.z = z;
	id = identity;

	Temperature = NULL;
}

void temperatureNode::allocateTemperatureTimeSeries(Real* Ttseries, long nentries)
{
	if (Temperature) return;

	Temperature = (Real*)calloc(nentries, sizeof(Real));

	if (Temperature)	memcpy(Temperature, &Ttseries[this->id*nentries], nentries * sizeof(Real));
}

void tempCellClass::allocateTemperatureTimeLocalNodes(int nentries)
{

	for (short i = c000; i <= c111; i++)
	{
		this->data.vertexTemp[i]->allocateTemperatureTimeSeries(owner()->nodes->temperatureSeriesRawData, nentries);
	}
}

tempCellClass::tempCellClass(tempCellPoolP own, TNodeP T000, TNodeP T001, TNodeP T010, TNodeP T011, TNodeP T100, TNodeP T101, TNodeP T110, TNodeP T111)
{
	initialCellNumber = 0;
	mySpace.xleft = 0;
	mySpace.xright = 0;
	mySpace.yleft = 0;
	mySpace.yright = 0;
	mySpace.zleft = 0;
	mySpace.zright = 0;
	data.vertexTemp[c000] = T000;
	data.vertexTemp[c001] = T001;
	data.vertexTemp[c010] = T010;
	data.vertexTemp[c011] = T011;
	data.vertexTemp[c100] = T100;
	data.vertexTemp[c101] = T101;
	data.vertexTemp[c110] = T110;
	data.vertexTemp[c111] = T111;
	this->own = own;
	next = NULL;
	prev = NULL;
}

tempCellClass::tempCellClass(tempCellPoolP owner, short xmin, short xmax, short ymin, short ymax, short zmin, short zmax)
{
	initialCellNumber = (zmax-zmin)*(ymax-ymin)*(xmax-xmin);
	mySpace.xleft = xmin;
	mySpace.xright = xmax;
	mySpace.yleft = ymin;
	mySpace.yright = ymax;
	mySpace.zleft = zmin;
	mySpace.zright = zmax;
	data.vertexTemp[c000] = NULL;
	data.vertexTemp[c001] = NULL;
	data.vertexTemp[c010] = NULL;
	data.vertexTemp[c011] = NULL;
	data.vertexTemp[c100] = NULL;
	data.vertexTemp[c101] = NULL;
	data.vertexTemp[c110] = NULL;
	data.vertexTemp[c111] = NULL;
	this->own = owner;

	if (own->first)	own->first->prev = this;
	this->next = own->first;
	this->prev = NULL;
	own->first = this;
	if (!(own->last))	own->last = this;
	id = 0;
}

tempCellClass::tempCellClass(void)
{
	initialCellNumber = 0;
	mySpace.xleft = 0;
	mySpace.xright = 0;
	mySpace.yleft = 0;
	mySpace.yright = 0;
	mySpace.zleft = 0;
	mySpace.zright = 0;
	data.vertexTemp[c000] = NULL;
	data.vertexTemp[c001] = NULL;
	data.vertexTemp[c010] = NULL;
	data.vertexTemp[c011] = NULL;
	data.vertexTemp[c100] = NULL;
	data.vertexTemp[c101] = NULL;
	data.vertexTemp[c110] = NULL;
	data.vertexTemp[c111] = NULL;
	this->own = NULL;
	next = NULL;
	prev = NULL;
}

Real tempCellClass::interpolateTemperatureSpace(short ix, short iy, short iz, Real currentTime, long timeKey, long nextKey)
{
	short x0 = mySpace.xleft;
	short xf = mySpace.xright;
	short y0 = mySpace.yleft;
	short yf = mySpace.yright;
	short z0 = mySpace.zleft;
	short zf = mySpace.zright;

	//cout << x0 << "\t" << y0 << "\t" << z0 << "\t" << xf << "\t" << yf << "\t" << zf << "\t" << endl; // for debug remove
	
	TNodeP TCellNodes[8] = { data.vertexTemp[c000],data.vertexTemp[c001],data.vertexTemp[c010],data.vertexTemp[c011],
		data.vertexTemp[c100],data.vertexTemp[c101],data.vertexTemp[c110],data.vertexTemp[c111] };

	Real currentNodesTemperature[8] = { 0 };

	Real xd = (Real)(ix - x0) / (Real)(xf - x0);
	Real yd = (Real)(iy - y0) / (Real)(yf - y0);
	Real zd = (Real)(iz - z0) / (Real)(zf - z0);

	//long timeKey = own->T_timer.key;
	
	 QUICKASSERT(timeKey >= 0 && timeKey < own->T_timer.size); //#as: it checks the condition. Must be uncommented after debugging
	 
	if (!(timeKey >= 0 && timeKey < own->T_timer.size)) //#as for debug
	//	cout << "here is the error" << endl;

	Real prevTime = 0;
	Real prevTime = own->T_timer.timeSeries[timeKey];

	Real nextTime;

	if (nextKey > timeKey)
	{
		if (timeKey < own->T_timer.size - 1) nextTime = own->T_timer.timeSeries[nextKey];
		else nextTime = own->T_timer.timeSeries[timeKey];
	}
	else if (nextKey < timeKey)
	{
		if (timeKey > 0) nextTime = own->T_timer.timeSeries[nextKey];
		else nextTime = own->T_timer.timeSeries[timeKey];
	}
	else
	{
		nextTime = own->T_timer.timeSeries[nextKey];
	}

//cout << prevTime << "\t" << currentTime << "\t" << nextTime << endl; //as for debug delete
	//QUICKASSERT(currentTime >= MINIMUM(prevTime, nextTime, nextTime) && currentTime <= MAXIMUM(prevTime, nextTime, nextTime));
	if (currentTime < MINIMUM(prevTime, nextTime, nextTime))
	{
		cout << prevTime << "min\t" << currentTime << "\t" << nextTime << endl; //as for debug delete
	}
	if (currentTime > MAXIMUM(prevTime, nextTime, nextTime))
	{
		cout << prevTime << "\tmax\t" << currentTime << "\t" << nextTime << endl; //as for debug delete
		
	}
	QUICKASSERT(currentTime >= MINIMUM(prevTime, nextTime, nextTime) && currentTime <= MAXIMUM(prevTime, nextTime, nextTime));
	//prevTime must correspond to time[timeKey] and nextTime to time[nextKey]; currentTime is the point between prev and next.

	for (int i = c000; i <= c111; i++)
		currentNodesTemperature[i] = TCellNodes[i]->interpolateTemperatureTime(prevTime, currentTime, nextTime, timeKey, nextKey);

	Real T00 = currentNodesTemperature[c000] * (1 - xd) + currentNodesTemperature[c100] * xd;
	Real T01 = currentNodesTemperature[c001] * (1 - xd) + currentNodesTemperature[c101] * xd;
	Real T10 = currentNodesTemperature[c010] * (1 - xd) + currentNodesTemperature[c110] * xd;
	Real T11 = currentNodesTemperature[c011] * (1 - xd) + currentNodesTemperature[c111] * xd;
	Real T0 = T00 * (1 - yd) + T10 * yd;
	Real T1 = T01 * (1 - yd) + T11 * yd;
	Real T = T0 * (1 - zd) + T1 * zd;

	return T;
}


// #as: modified function
Real tempCellClass::determineNucleationTime(short ix, short iy, short iz)
{
#define REGDATSIZ (5)
	//#define MELTINGPOINT (1375) //#as reading from input file
	int check = 0;
	Real Tcrit = 0;
	/*Real Tt = 0;*/
	long ti = 0;
	long Tem = 0;
	Real Tmax = 0;
	caHdlP ca = automaton();
	tempCellPoolP TCP = owner();
	Real* time = TCP->T_timer.timeSeries;
	long len = TCP->T_timer.size;
	Real MELTINGPOINT = ca->gen.Tmelt;
	Real t[5] = { 0 };
	Real T[5] = { 0 };
	Real m = 0; //slope
	Real c0 = 0; //intercept
	Real R2 = 0; // R-squared
	//Real mm = 0; //slope
	//Real cc0 = 0; //intercept
	//Real RR2 = 0; // R-squared
	Real meltTime = NO_MELT;
	TCP->shortesttNucleationtime = 9999; 
	/*Real tt[] = { 0 };
	Real TT[] = { 0 };*/
	//caHdl h;
	for (long i = len - 1; i >= 0; i--)
	{
		long nextPos = 0;
		if (i == 0)
			nextPos = 0;

		else nextPos = i - 1;
		Real Tem = interpolateTemperatureSpace(ix, iy, iz, time[i], i, nextPos);
		if (check == 0) {
			if (Tem > MELTINGPOINT)
			{
				Tcrit = Tem;
				check = 1;
				ti = i;

			}
		}
	
		if (Tem >= Tmax) {
			Tmax = Tem;
		}
		
	}
	if (Tmax > 2 * MELTINGPOINT)
		return NO_MELT; //overmelted

	if (Tmax < MELTINGPOINT)
	{
		/*long index = (iz - 0) * ca->cells.data.zFac + (iy - 0) * ca->cells.data.xPer + (ix - 0);
		a.total[index]= (Real)0x00000001;*/
		ca->PT.tmax.push_back(Tmax);
		ca->PT.x.push_back(ix);
		ca->PT.y.push_back(iy);
		ca->PT.z.push_back(iz);

		//maxtemp.tmax = Tmax;
		//h.nomelt.push_back(maxtemp);
		//cout << "Tmax is " << Tmax << endl;
		return NO_MELT;
	}
	for (long j = 0; j < REGDATSIZ; j++)
	{
		long nextPos = 0;
		//long nextPosp = 0;
		if (ti - j == 0) nextPos = 0;  
		else nextPos =ti- j - 1;
		//if (ti + j == len-1) nextPosp = 0; 
			//else nextPosp = j + 1;
		t[j] = own->T_timer.timeSeries[ti - j];
		//tt[j] = own->T_timer.timeSeries[ti + j];
		T[j] = interpolateTemperatureSpace(ix, iy, iz, t[j], ti - j, nextPos);

		//TT[j] = interpolateTemperatureSpace(ix, iy, iz, t[j], ti + j, nextPos);
		
		//cout << "coolinng temp is:" << T[j] <<"at time = "<< t[j] << endl;
	}
	if (ti == len - 1 && T[0] > MELTINGPOINT) // if it remains liquid
		return LIQUID;
	else
		//ca->linearRegression(tt, TT, REGDATSIZ, &mm, &cc0, &RR2);
		ca->linearRegression(t, T, REGDATSIZ, &m, &c0, &R2);
	//cout << "m :" << m << "c:" << c0 << endl;//#as: for debug
	//cout << "mm :" << mm << "cc:" << cc0 << endl;//#as: for debug

	//Real Tcrit = m * t[REGDATSIZ - 1] + c0;
	//Real meltTime1 = (Tmax - c0) / m;
	//cout << "Tmax melt time is " << meltTime1 << endl;
	//Real meltTime2 = (MELTINGPOINT - cc0) / mm;
	//cout << "meltTime2 is " << meltTime << endl;
	//cout << "Tcrit is" << Tcrit << endl;
	//meltTime = (Tcrit - c0) / m;
	//cout << "melt time is " << meltTime << endl;
	meltTime = (MELTINGPOINT - c0) / m;
	if (meltTime < TCP->shortesttNucleationtime)
		TCP->shortesttNucleationtime = meltTime;
	//cout << "melt time isas " << meltTime << endl;
	return meltTime;


}


// originally 
/*
Real tempCellClass::determineNucleationTime(short ix, short iy, short iz)
{
#define REGDATSIZ (5)
//#define MELTINGPOINT (1375) //#as reading from input file

	caHdlP ca = automaton();
	tempCellPoolP TCP = owner();
	Real* time = TCP->T_timer.timeSeries;
	long len = TCP->T_timer.size;
	Real MELTINGPOINT = ca->gen.Tmelt;
	Real t[5] = { 0 };
	Real T[5] = { 0 };
	Real m = 0; //slope
	Real c0 = 0; //intercept
	Real R2 = 0; // R-squared
	Real meltTime = NO_MELT;

	for (long i = len - 1; i > REGDATSIZ; i--)
	{
	//	cout << "i:" << i << endl; //#as debug
	//	cout << "t:" << (*time) << endl;//#as debug
		for (long j = 0; j < REGDATSIZ; j++)
		{
			long nextPos = 0;

			if (i - j == 0) nextPos = 0;
			else nextPos = i - j - 1;

			t[j] = time[i - j];
			T[j] = interpolateTemperatureSpace(ix, iy, iz, t[j], i-j, nextPos);
			// cout << i << "\t\t" << j << "\t\t" << "T[j]" << T[j] << "\t" << "time:" << t[j] << endl; #as debug
		}

		if (i == len - 1 && T[0] > MELTINGPOINT && T[REGDATSIZ] > MELTINGPOINT)
			return LIQUID;

		ca->linearRegression(t, T, REGDATSIZ, &m, &c0, &R2);

		Real Tcrit = m * t[REGDATSIZ - 1] + c0;
		//cout << "Tcrit: "<<Tcrit << endl; //debug
		if (Tcrit <= MELTINGPOINT) 
			continue;
		else
		{
			meltTime = (MELTINGPOINT - c0) / m;
			return meltTime;
		}

	}
	return NO_MELT ; // change to liquid
	//return meltTime;
}*/

Real temperatureNode::interpolateTemperatureTime(Real t0, Real tc, Real tf, long pos, long nextPos)
{
	Real T0 = Temperature[pos];
	Real Tf = Temperature[nextPos];
	
	//if (tf == t0 == tc && T0 == Tf) 
	//{
	//	return Tf;
	//}

	Real Tc = T0 + (Tf - T0) * ((tc - t0) / (tf - t0));
	return Tc;
}

tempNodePool::tempNodePool(tempCellPool* own)
{
	this->owner = own;
}

tempNodePool::~tempNodePool(void)
{
	for (int i = 0; i < this->size; i++) delete this->all[i];
	if (all) free(all);
}

//for reading txt file


//void tempNodePool::readNodalTemperatureTimeSeries(const char* filename)
//{
//	long rows = owner->T_timer.size;
//	long columns = this->size;
//
//	temperatureSeriesRawData = (Real*)calloc((size_t)rows * (size_t)columns, sizeof(Real));
//
//	if (!temperatureSeriesRawData) owner->owner->err->reportError(ERRTXT("Cannot allocate more memory"));
//
//	nodalTemperatureFilename = new char[BUFSIZ];
//
//	strcpy(nodalTemperatureFilename, filename);
//
//	ifstream input(nodalTemperatureFilename);
//	string line;
//	int l = 0;
//	long counter = 0;
//
//	while (getline(input, line))
//	{
//		float value;
//
//		stringstream ss(line);
//
//		int k = 0;
//
//		//ss.seekg(0, ios::end);
//		//long size = ss.tellg();
//
//		while (ss >> value)
//		{
//			temperatureSeriesRawData[k+k*(rows-1)+l] = value;
//			++k;
//			++counter;
//		}
//		++l;
//	}
//	input.close();
//	delete nodalTemperatureFilename;
//	nodalTemperatureFilename=NULL;
//	QUICKASSERT(counter==rows*columns);
//}




// for reading binary file


void tempNodePool::readNodalTemperatureTimeSeries(const char* filename)
{
	long rows = owner->T_timer.size;
	long columns = this->size;

	temperatureSeriesRawData = (Real*)calloc((size_t)rows * (size_t)columns, sizeof(Real));

	if (!temperatureSeriesRawData) owner->owner->err->reportError(ERRTXT("Cannot allocate more memory"));

	nodalTemperatureFilename = new char[BUFSIZ];

	strcpy(nodalTemperatureFilename, filename);

	ifstream input(nodalTemperatureFilename, ios::in |ios::binary);

	uint32_t integer_value;
	long counter = 0;
	int l = 0;
	input.read(reinterpret_cast<char*>(&integer_value), sizeof(integer_value));

	while (!input.eof())
	{
		int k = 0; 
		while (k < columns)
		{
			temperatureSeriesRawData[k + k * (rows - 1) + l] = Real(integer_value) / 1000;
			input.read(reinterpret_cast<char*>(&integer_value), sizeof(integer_value));
			++k;
			++counter;
		}
		++l;
	}

	input.close();
	delete nodalTemperatureFilename;
	nodalTemperatureFilename = NULL;
	QUICKASSERT(counter == rows * columns);
}


tempCellClass::~ tempCellClass(void)
{

}







Real tempCellClass::TempParticle(short ix, short iy, short iz) //defined in line 117 in tempcontrol.h
{
	#define REGDATSIZ (5)
	//#define MELTINGPOINT (1375) //#as reading from input file
	int check = 0;
	Real Tcrit = 0;
	/*Real Tt = 0;*/
	long ti = 0;
	long Tem = 0;
	Real Tmax = 0;
	caHdlP ca = automaton();
	tempCellPoolP TCP = owner();
	Real* time = TCP->T_timer.timeSeries;
	long len = TCP->T_timer.size;
	Real MELTINGPOINT = ca->gen.Tmelt;
	Real t[5] = { 0 };
	Real T[5] = { 0 };
	Real m = 0; //slope
	Real c0 = 0; //intercept
	Real R2 = 0; // R-squared
	
	Real meltTime = NO_MELT;
	
	for (long i = len - 1; i >= 0; i--)
	{
		long nextPos = 0;
		if (i == 0)
			nextPos = 0;

		else nextPos = i - 1;
		Real Tem = interpolateTemperatureSpace(ix, iy, iz, time[i], i, nextPos);
		if (check == 0) {
			if (Tem > MELTINGPOINT)
			{
				Tcrit = Tem;
				check = 1;
				ti = i;

			}
		}
		
		if (Tem >= Tmax)
		{
			Tmax = Tem;
		}
		//return Tcrit;//#as: for debug
	}
	return Tmax;
}


	/*
Real tempCellClass::TempParticle(short xl, short xr, short yl, short yr, short zl, short zr, short ix, short iy, short iz, Real currentTime, long timeKey, long nextKey)
// this fucntion is declared in tempcontrol.h line 161
{
	short x0 = xl;
	short xf = xr;
	short y0 = yl;
	short yf =yr;
	short z0 = zl;
	short zf = zr;

	//cout << x0 << "\t" << y0 << "\t" << z0 << "\t" << xf << "\t" << yf << "\t" << zf << "\t" << endl; // for debug remove

	TNodeP TCellNodes[8] = { data.vertexTemp[c000],data.vertexTemp[c001],data.vertexTemp[c010],data.vertexTemp[c011],
		data.vertexTemp[c100],data.vertexTemp[c101],data.vertexTemp[c110],data.vertexTemp[c111] };

	Real currentNodesTemperature[8] = { 0 };

	Real xd = (Real)(ix - x0) / (Real)(xf - x0);
	Real yd = (Real)(iy - y0) / (Real)(yf - y0);
	Real zd = (Real)(iz - z0) / (Real)(zf - z0);

	//long timeKey = own->T_timer.key;

	QUICKASSERT(timeKey >= 0 && timeKey < own->T_timer.size);



	Real prevTime = 0;
	prevTime = own->T_timer.timeSeries[timeKey];

	Real nextTime;

	if (nextKey > timeKey)
	{
		if (timeKey < own->T_timer.size - 1) nextTime = own->T_timer.timeSeries[nextKey];
		else nextTime = own->T_timer.timeSeries[timeKey];
	}
	else if (nextKey < timeKey)
	{
		if (timeKey > 0) nextTime = own->T_timer.timeSeries[nextKey];
		else nextTime = own->T_timer.timeSeries[timeKey];
	}
	else
	{
		nextTime = own->T_timer.timeSeries[nextKey];
	}

	//cout << prevTime << "\t" << currentTime << "\t" << nextTime << endl; //as for debug delete
	QUICKASSERT(currentTime >= MINIMUM(prevTime, nextTime, nextTime) && currentTime <= MAXIMUM(prevTime, nextTime, nextTime));

	//prevTime must correspond to time[timeKey] and nextTime to time[nextKey]; currentTime is the point between prev and next.

	for (int i = c000; i <= c111; i++)
		currentNodesTemperature[i] = TCellNodes[i]->interpolateTemperatureTime(prevTime, currentTime, nextTime, timeKey, nextKey);

	Real T00 = currentNodesTemperature[c000] * (1 - xd) + currentNodesTemperature[c100] * xd;
	Real T01 = currentNodesTemperature[c001] * (1 - xd) + currentNodesTemperature[c101] * xd;
	Real T10 = currentNodesTemperature[c010] * (1 - xd) + currentNodesTemperature[c110] * xd;
	Real T11 = currentNodesTemperature[c011] * (1 - xd) + currentNodesTemperature[c111] * xd;
	Real T0 = T00 * (1 - yd) + T10 * yd;
	Real T1 = T01 * (1 - yd) + T11 * yd;
	Real T = T0 * (1 - zd) + T1 * zd;

	return T;
}
*/
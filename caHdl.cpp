//#define _CRT_SECURE_NO_WARNINGS //as i hvae disbled pole figure line 4020 pole figure LINE 2296 and 2249 2306

#include <string.h>
#include <time.h>
#include "caHdl.h"
#include <math.h>
#include <ctime>
#include <stdio.h>
#include <time.h>
#include "io.h"

caHdl::~caHdl(void) //individual cells need to be deleted ??
{
	delete parallelHdl;
	delete err;

	fclose(core8);
	free(this->cells.all);
	if (theFile) MPI_File_close(&theFile);
	delete Tcells;
	endParallelProcesses();
	//seeds in destructor?

	ostringstream message;
	std::time_t time1 = std::time(nullptr);// #added rg
	message << "Parallel Processes Ends" << std::asctime(std::localtime(&time1));// #added rg
	write(message.str(), &logFile);// #added rg

	logFile.close();

}

caHdl::caHdl(void) : polyCrystal(this)
{
	myRandom.init(-46356);//-46356 : original seed value.
	//myRandom.init (-currentTime(NULL));
	this->time = 0.0;
	this->deltaTime = 0.0;
	grains.owner = this;
	//	boundaries.owner = this;
	err = new error;

	// tempCounter = 0; // No use.

	Real constants[27] = { FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE,
					  FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE, FGEOFACE, FGEOEDGE, FGEOFACE,
					  FGEOEDGE,	FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG,
					  FGEOEDGE, FGEODIAG, 1 };
					/* {FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE,
					  FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEOFACE, FGEOFACE, FGEOEDGE, FGEOFACE,
					  FGEOEDGE,	FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG,
					  FGEOEDGE, FGEODIAG, 1 };*/
					/* {FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG, FGEOEDGE,
					  FGEODIAG, .707, .5, .707, 1, 1, .707, .5,.707,
					  FGEODIAG, FGEOEDGE, FGEODIAG, FGEOEDGE, FGEOFACE, FGEOEDGE, FGEODIAG,
					  FGEOEDGE, FGEODIAG, 1 };*/
		

	memcpy(parm.fgeo, constants, 27 * sizeof(Real));

	Tcells = new tempCellPool(this);
	theFile = NULL;
}

void caHdl::initializeSolCA(const char* inputfile, int iautomaton)
{
	Real Ncells;

	printf("Initializing...\n");

	//err = new error;

	ostringstream logCAfname, logfname;
	logCAfname << "logCA.";
	logfname << "log.";

	//linux without leading zeros

	logCAfname << iautomaton << ".dat";
	logfname << iautomaton << ".csv";

	if (!open(logCAfname.str().c_str(), &logFile)) err->reportError(ERRTXT("Cannot open logFile"));

	ostringstream message;// #added rg
	std::time_t time1 = std::time(nullptr);// #added rg
	message << "Start time is " << std::asctime(std::localtime(&time1));// #added rg
	write(message.str(), &logFile);// #added rg

	err->logFile = &logFile;

	if (!open(READ, inputfile, &(this->core8)))
	{
		err->reportError(ERRTXT("Cannot open Input. Termination..."));

		endProgram();												///endProgram() Not yet implemented
	}

	getGeneralInput(&Ncells);

	determineCAsizeS(Ncells);

	initializeGrainsSolCA(); // Grain initialization. Here we assign entire domain as one grain of orientation 0,0,0.
}

void caHdl::getGeneralInput(Real* Ncells)
{

	dataBlockP temperatureTimeData; // reads the time steps.
	temperatureTimeData = readDataBlock("TemperatureTime", core8); // what is core 8
	readTemperatureTimeDiscretization(temperatureTimeData);
	delete temperatureTimeData;

	dataBlockP temperatureNodes; //reads the temp node coordiantes.
	temperatureNodes = readDataBlock("TemperatureNodeCoordinates", core8);
	this->Tcells->nodes->tnodesraw = readTemperatureNodes(temperatureNodes);
	delete temperatureNodes;

	dataBlockP temperatureCells; // reads the temp_cell neighbouring temp_node id
	temperatureCells = readDataBlock("TemperatureCells", core8);
	this->Tcells->tcellraw = readTemperatureCells(temperatureCells);
	delete temperatureCells;

	dataBlockP TC_x_discr; // x coordiantes of temp_cells
	dataBlockP TC_y_discr;
	dataBlockP TC_z_discr;
	TC_x_discr = readDataBlock("TemperatureCellsDiscretizationX", core8);
	TC_y_discr = readDataBlock("TemperatureCellsDiscretizationY", core8);
	TC_z_discr = readDataBlock("TemperatureCellsDiscretizationZ", core8);
	readTemperatureCellsDiscretization(TC_x_discr, TC_y_discr, TC_z_discr);
	delete TC_x_discr;
	delete TC_y_discr;
	delete TC_z_discr;

	dataBlockP parmDataNucGrains;
	parmDataNucGrains = readDataBlock("ListOfNucleusOrientations", core8);
	numberInputNucOris = parmDataNucGrains->lineCount;
	grainOriNuc = (double*)calloc(numberInputNucOris * 3, sizeof(double));
	dataLineP lineNucGrains = parmDataNucGrains->first;
	for (int ing = 0; ing < numberInputNucOris; ing++)
	{
		int idx = ing * 3; //why this multiplication also in this is present in bulk nucleation
		grainOriNuc[idx] = _PI_ * getReal(lineNucGrains, 1) / 180.;
		grainOriNuc[idx + 1] = _PI_ * getReal(lineNucGrains, 2) / 180.;
		grainOriNuc[idx + 2] = _PI_ * getReal(lineNucGrains, 3) / 180.;

		lineNucGrains = lineNucGrains->next;
	}
	delete lineNucGrains;
	delete parmDataNucGrains;



	dataBlockP parmData;
	parmData = readDataBlock("ParameterList", core8);

	Tcells->nodes->readNodalTemperatureTimeSeries(geTString("NodalTemperatureTimeSeriesFilename", parmData));
	*Ncells = geTReal("NumberOfCaCellsInMio", parmData) * 1.0e6; //#as:number of total cells 
	gen.Uc = geTReal("criticalundercooling", parmData);   //#as check for undercooling from input file
	gen.Tmelt = geTReal("MeltingTemperature", parmData);//#as melting temperarure  from input file
	gen.eta1 = geTReal("eta1", parmData);//#as melting temperarure  from input file
	gen.eta2 = geTReal("eta2", parmData);//#as melting temperarure  from input file
	gen.Tcrit_mean= geTReal("CriticalMeanTemp", parmData);//#as melting temperarure  from input file
	gen.Tcrit_sigma= geTReal("CriticalTempDeviation", parmData);//#as melting temperarure  from input file





	gen.ASPIREDFILLPERSTEP = ASPIREDFILLPERSTEP_DEFAULT;//set DEFAULT ASPIREDFILLPERSTEP and MAXFILLPERSTEP
	gen.MAXFILLPERSTEP = MAXFILLPERSTEP_DEFAULT; //is at least initialized, values used in Carmen Schaefers work, good compromise between speed and accuracy

	gen.ASPIREDFILLPERSTEP = geTReal("ASPIREDFILLPERSTEP", parmData);
	gen.MAXFILLPERSTEP = geTReal("MAXFILLPERSTEP", parmData);

	if (gen.ASPIREDFILLPERSTEP > 1 || gen.ASPIREDFILLPERSTEP <= 0)
		err->reportError(ERRTXT("ASPIREDFILLPERSTEP needs to be in the interval 0 < ASPIREDFILLPERSTEP <= 1."));

	if (gen.MAXFILLPERSTEP > 1 || gen.MAXFILLPERSTEP <= 0)
		err->reportError(ERRTXT("MAXFILLPERSTEP needs to be in the interval 0 < MAXFILLPERSTEP <= 1."));

	if (gen.ASPIREDFILLPERSTEP > gen.MAXFILLPERSTEP)
		err->reportError(ERRTXT("Severe mistakes are expected if the ASPIREDFILLSTEP is large than the MAXFILLPERSTEP."));

	//protocols of periodicity and communication for parallel computations
	char* periodicity = geTString("Periodicity", parmData);
	if (stringsEqual(periodicity, "boxPer"))				this->boxPer = BOXPER;
	else if (stringsEqual(periodicity, "nonBoxPer"))	this->boxPer = NONBOXPER;
	else													this->boxPer = BOXPER;

	char* nodePeriodicity = geTString("nodePeriodicity", parmData);
	if (stringsEqual(nodePeriodicity, "nodePer"))		this->nodePer = NODE_PER;
	else if (stringsEqual(nodePeriodicity, "nonNodePer"))	this->nodePer = NON_NODE_PER;
	else													this->nodePer = NON_NODE_PER;

	char* protocolType = geTString("CommunicationProtocol", parmData);
	if (stringsEqual(protocolType, "noComm"))				this->protocol = NOCOMM;
	else if (stringsEqual(protocolType, "deterministic"))	this->protocol = DETERMINISTIC;
	else if (stringsEqual(protocolType, "statistical"))		this->protocol = STATISTICAL;
	else	err->reportError(ERRTXT("Communication protocol not recognized. Valid options: noComm, deterministic and statistical"));

	if (protocol == STATISTICAL)
	{
		char* transferType = geTString("TransferType", parmData);
		transfer.frequency = (long)geTReal("TransferFrequency", parmData);

		if (stringsEqual(transferType, "steps"))			this->commType = STEPS;
		else	err->reportError(ERRTXT("comm: statistical but no transfer type recognized. Valid options: steps"));

		char* insertion = geTString("NucleiInsertionSite", parmData); // here is problem. "NucleiInsertionSite" is not defined in input file.

		if (stringsEqual(insertion, "bulk"))				this->insertionSite = BULK;
		else if (stringsEqual(insertion, "interface"))		this->insertionSite = INTERFACE;
		else	err->reportError(ERRTXT("Insertion site for nuclei not recognized. Valid options: bulk and interface"));
	}

	if (!(COMM_ENABLED & protocol))
		transfer.never = NEVER;

	// Grain boundary character and intrinsic mobilities
	parm.TwistAngleDeviation = 5.0; // Not needed
	parm.AxisAngleDeviation = 5.0; // Not needed
	parm.HighMobAngle = 40.0;	//Growth Selection 40∞<111>

	// Output parameters.
	mstOutput.frequency = geTInt("PlotFrequencyMicroStructure", parmData);
	txtOutput.frequency = geTInt("PlotFrequencyTexture", parmData);
	rxOutput.frequency = geTInt("PlotFrequencyRxFrac", parmData);
	rvOutput.frequency = geTInt("PlotFrequencyRVoutput", parmData);

	//determine desired user output strategy
	char* msoutput = geTString("MSOutput", parmData);
	if (strcmp(msoutput, "user") == 0)	userDefinedMSOutput = true;
	else								userDefinedMSOutput = false;

	char* msoutputbinaryopt = geTString("GenerateMSBinaryOption", parmData);
	if (strcmp(msoutputbinaryopt, "yes") == 0) {
		generateMSBinary = true;

		char* msplotmode = geTString("PlotMode", parmData);
		if (stringsEqual(msplotmode, "inf"))		this->MSBinaryPlotMode = PLOTINFPROGRESS;
		else if (stringsEqual(msplotmode, "ori"))	this->MSBinaryPlotMode = PLOTRGB;
		else										err->reportError(ERRTXT("MSBinary Output desired but no valid option specified under parmData parameter PlotMode - valid options: inf, ori"));
	}
	else generateMSBinary = false;

	char* backgroundMode = geTString("Background", parmData);

	if (stringsEqual(backgroundMode, "white"))			gen.backgroundMode = WHITE;
	else if (stringsEqual(backgroundMode, "black"))		gen.backgroundMode = BLACK;
	else if (stringsEqual(backgroundMode, "deformed"))	gen.backgroundMode = DEFMODUS;
	else													err->reportError(ERRTXT("Plotting background not recognized"));

	if (userDefinedMSOutput)
	{
		dataBlockP userDefinedFData;
		userDefinedFData = readDataBlock("UserDefinedMSOutputFrequency", core8);
		dataLineP line = userDefinedFData->first;
		while (line) {
			char* parameter = getString(line, 1);
			if (strcmp("\"X\"", parameter) == 0) {
				Real value = getReal(line, 2);
				if (value >= 0 && value <= 1) {
					userDefinedFrequencies.push_back(value);
				}
				else err->reportError(ERRTXT("User defined Frequencies must be between 0 and 1"));
			}
			line = line->next;
		}
		userDefinedFrequencies.sort(); // sort list by frequency
		delete line;
		delete userDefinedFData;
	}
	outCount = 0;

	dataBlockP userDefinedOriFData;
	userDefinedOriFData = readDataBlock("UserDefinedOriOutputFrequency", core8);
	dataLineP oriLine = userDefinedOriFData->first;
	while (oriLine) {
		char* parameter = getString(oriLine, 1);
		if (strcmp("\"X\"", parameter) == 0) {
			Real value = getReal(oriLine, 2);
			if (value >= 0 && value <= 1) {
				userDefinedOriFrequencies.push_back(value);
			}
			else err->reportError(ERRTXT("User defined Frequencies must be between 0 and 1"));
		}
		oriLine = oriLine->next;
	}
	userDefinedOriFrequencies.sort();
	oriOutCount = 0;
	delete oriLine;
	delete userDefinedOriFData;

	delete parmData;

	ostringstream message;
	message << "ASPIREDFILLPERSTEP:" << gen.ASPIREDFILLPERSTEP << endl;
	message << "MAXFILLPERSTEP:" << gen.MAXFILLPERSTEP << endl;
	message << "TIMESLOTCOUNT:" << TIMESLOTCOUNT << endl;
	message << "FGEOFACE:" << FGEOFACE << endl;
	message << "FGEOEDGE:" << FGEOEDGE << endl;
	message << "FGEODIAG:" << FGEODIAG << endl;

	message << "Melting Temperature:" << gen.Tmelt << endl; //#as reading from input file
	message << "Critical undercooling:" << gen.Uc << endl;//#as reading from input file
	message << "eta 1: " << gen.eta1 << endl; //#as reading from input file
	message << "eta 2: " << gen.eta2 << endl; //#as reading from input file

	message << "Texture Output:" << txtOutput.frequency << " steps" << endl;
	message << "Rx Output:" << rxOutput.frequency << " steps" << endl;
	message << "Microstructure Output:" << mstOutput.frequency << " steps" << endl;
	message << "Recovery Output:" << rvOutput.frequency << " steps" << endl << endl; //Not needed.
	message << "HighMobAngle:" << parm.HighMobAngle << endl; // Not needed
	message << "TwistAngleDeviation:" << parm.TwistAngleDeviation << endl; // Not neeeded
	message << "AxisAngleDeviation:" << parm.AxisAngleDeviation << endl; // Not needed
	write(message.str(), &logFile);
}

void caHdl::readTemperatureTimeDiscretization(dataBlockP tTData) // time points assignments
{

	long Nentries = tTData->lineCount;
	Tcells->T_timer.key = 0;
	Tcells->T_timer.size = Nentries;
	Tcells->T_timer.timeSeries = (Real*)calloc(Nentries, sizeof(Real));

	dataLineP lineTimeTemperature = tTData->first; // memory 

	double dlt = 5;

	//ofstream check;//**added to check the output
	//check.open("readTemperatureTimeDiscretization.txt");//**added to check the output
	//check << "Nentries are:" << Nentries << "\n \n"<< endl; //**added to check the output
	//check << "lineTimeTemperature is:" << (lineTimeTemperature) <<"\n\n"<< endl;//**added to check the output
	//check << "i\t\t\t"<<"lineTimeTemperatureseries(i)\t\t\t\t\t\t lineTimeTemperature" << endl; //**added to check the output


	for (int i = 0; i < Nentries; i++)
	{

		Tcells->T_timer.timeSeries[i] = getReal(lineTimeTemperature, 1);
		lineTimeTemperature = lineTimeTemperature->next;

		//*****
		//check << i << "\t\t\t" << Tcells->T_timer.timeSeries[i]<<"\t\t\t\t\t\t"<< lineTimeTemperature << endl;//**added to check the output


	}

	Tcells->setTimeEncodeFactor();
}

TNodeRawData* caHdl::readTemperatureNodes(dataBlockP Tnodes)
{


	long Nentries = Tnodes->lineCount;
	TNodeRawData* nodesData = (TNodeRawData*)calloc(Nentries, sizeof(TNodeRawData));

	dataLineP lineTNodes = Tnodes->first;

	Real xmin = 1e10;
	Real xmax = 0.;
	Real ymin = 1e10;
	Real ymax = 0.;
	Real zmin = 1e10;
	Real zmax = 0.;

	//************added to see ouput
	//ofstream check;//************added to see ouput
	//check.open("readTemperatureNodes.txt");//************added to see ouput
	//check << "valye of Nentries are:" << Nentries << endl;//************added to see ouput
	//check << "(i)\t\t\t\t\t\t id\t\t\t\t\t\t posx\t\t\t\t\t\t posy\t\t\t\t\t\t posz\t\t\t\t\t\t LineTnodes"<<endl;//************added to see ouput
	//************added to see ouput


	for (int i = 0; i < Nentries; i++)
	{
		long id = getInt(lineTNodes, 1) - 1;

		nodesData[i].id = id;

		Real x = getReal(lineTNodes, 2);
		Real y = getReal(lineTNodes, 3);
		Real z = getReal(lineTNodes, 4);

		if (x < xmin) xmin = x;
		if (x > xmax) xmax = x;
		if (y < ymin) ymin = y;
		if (y > ymax) ymax = y;
		if (z < zmin) zmin = z;
		if (z > zmax) zmax = z;

		nodesData[i].pos.x = x;
		nodesData[i].pos.y = y;
		nodesData[i].pos.z = z;

		lineTNodes = lineTNodes->next;

		//************added to see ouput

		//check << i << "\t\t\t\t\t\t" << nodesData[i].id<<"\t\t\t\t\t "<< nodesData[i].pos.x << "\t\t\t\t\t\t " << nodesData[i].pos.y << "\t\t\t\t\t\t " << nodesData[i].pos.z << "\t\t\t\t\t\t " << lineTNodes << endl; //added for output


		//************added to see ouput

	}

	this->Tcells->xmin = xmin;
	this->Tcells->xmax = xmax;
	this->Tcells->ymin = ymin;
	this->Tcells->ymax = ymax;
	this->Tcells->zmin = zmin;
	this->Tcells->zmax = zmax;

	this->Tcells->nodes->size = Nentries;

	///****dded to see ouput
	//check << "\n\n xmin\t ymin and \t zmin are:" << this->Tcells->xmin << "\t" << this->Tcells->ymin << "\t" << this->Tcells->zmin << "\t" << endl;
	//check.close();
	//*****dded to see ouput

	return nodesData;
}

void caHdl::readTemperatureCellsDiscretization(dataBlockP xdiscr, dataBlockP ydiscr, dataBlockP zdiscr)
{
	long nx = xdiscr->lineCount;
	long ny = ydiscr->lineCount;
	long nz = zdiscr->lineCount;

	this->Tcells->xPerT = (short)(nx - 1);
	this->Tcells->yPerT = (short)(ny - 1);
	this->Tcells->zPerT = (short)(nz - 1);

	this->Tcells->ntcells = this->Tcells->xPerT * this->Tcells->yPerT * this->Tcells->zPerT;

	this->Tcells->xdiscr = (Real*)calloc(nx, sizeof(Real));
	this->Tcells->ydiscr = (Real*)calloc(ny, sizeof(Real));
	this->Tcells->zdiscr = (Real*)calloc(nz, sizeof(Real));

	dataLineP lineDiscr = xdiscr->first;

	///****** added for output

	/*ofstream check;///****** added for output
	check.open("readTemperatureCellsDiscretization.txt");///****** added for output
	check << "values of nx\t ny\t and nz\t are:" << nx << "\t" << ny << "\t" << nz << endl;///****** added for output
	check << "values of xPerT\t yPerT\t and zPerT\t are:" << (nx-1) << "\t" << (ny-1) << "\t" << (nz-1) << endl;///****** added for output
	check << "(i)\t lineDiscr\t xdiscr\t" << endl;///****** added for output*/

	///****** added for output

	for (int i = 0; i < nx; i++)
	{
		this->Tcells->xdiscr[i] = getReal(lineDiscr, 1); //lineDiscr is name column is 1
		lineDiscr = lineDiscr->next;
		// check << i <<"\t\t"<< lineDiscr<<"\t\t" << this->Tcells->xdiscr[i] << "\t\t"<< endl;  ///****** added for output
	}

	//check << "(i)\t lineDiscr\t ydiscr\t" << endl;///****** added for output

	lineDiscr = ydiscr->first;



	for (int i = 0; i < ny; i++)
	{
		this->Tcells->ydiscr[i] = getReal(lineDiscr, 1);
		lineDiscr = lineDiscr->next;
		//check << i << "\t\t" << lineDiscr << "\t\t" << this->Tcells->ydiscr[i] << "\t\t" << endl;  ///****** added for output
	}
	//check << "(i)\t lineDiscr\t zdiscr\t" << endl;///****** added for output
	lineDiscr = zdiscr->first;

	for (int i = 0; i < nz; i++)
	{
		this->Tcells->zdiscr[i] = getReal(lineDiscr, 1);
		//check << i << "\t\t" << lineDiscr << "\t\t" << this->Tcells->zdiscr[i] << "\t\t" << endl;  ///****** added for output
		lineDiscr = lineDiscr->next;
	}
}

cellNodesRawData* caHdl::readTemperatureCells(dataBlockP Tcells)
{


	long Nentries = Tcells->lineCount;
	cellNodesRawData* tcells = (cellNodesRawData*)calloc(Nentries, sizeof(cellNodesRawData));

	dataLineP lineTcells = Tcells->first;

	/*///****** added for output
	ofstream check;///****** added for output
	check.open("readTemperatureCells.txt");///****** added for output
	check << "values of Nentries are:" << Nentries << endl;///****** added for output
	check << "i\t\t\t c000\t\t\t c001\t\t\t c010\t\t\t c011\t\t\t c100\t\t\t c101\t\t\t\ c110\t\t\t c111" << endl;
	///****** added for output*/


	for (int i = 0; i < Nentries; i++)
	{
		tcells[i].cellnodes[c000] = getInt(lineTcells, 2) - 1;
		tcells[i].cellnodes[c001] = getInt(lineTcells, 3) - 1;
		tcells[i].cellnodes[c010] = getInt(lineTcells, 4) - 1;
		tcells[i].cellnodes[c011] = getInt(lineTcells, 5) - 1;
		tcells[i].cellnodes[c100] = getInt(lineTcells, 6) - 1;
		tcells[i].cellnodes[c101] = getInt(lineTcells, 7) - 1;
		tcells[i].cellnodes[c110] = getInt(lineTcells, 8) - 1;
		tcells[i].cellnodes[c111] = getInt(lineTcells, 9) - 1;
		lineTcells = lineTcells->next;

		//check << i << "\t\t\t " << tcells[i].cellnodes[c000] << "\t\t\t " << tcells[i].cellnodes[c001] //as:added for output
		//		   << "\t\t\t " << tcells[i].cellnodes[c010] << "\t\t\t " << tcells[i].cellnodes[c011] //as:added for output
		//		   << "\t\t\t " << tcells[i].cellnodes[c100] << "\t\t\t " << tcells[i].cellnodes[c101]//as:added for output
		//		   << "\t\t\t " << tcells[i].cellnodes[c110] << "\t\t\t " << tcells[i].cellnodes[c111] << endl; //as:added for output

	}
	return tcells;
	//check.close();
}

void caHdl::determineCAsizeS(Real Ncells)
{
	ostringstream message;

	Real xmin = this->Tcells->xmin;
	Real xmax = this->Tcells->xmax;
	Real ymin = this->Tcells->ymin;
	Real ymax = this->Tcells->ymax;
	Real zmin = this->Tcells->zmin;
	Real zmax = this->Tcells->zmax;

	long cellsxyz = (long)(pow((double)Ncells, (double)1 / 3.0f) + 0.5); // number of cells in each direction from wherre it is coming
	long ncells = CUBE(cellsxyz); // total number of cells.

	Real Vtotal = (xmax - xmin) * (ymax - ymin) * (zmax - zmin); //..000000064
	Real lengthxyz = 0;
	lengthxyz = pow(Vtotal, (double)1 / 3.0f);//.004m

	this->xsize = this->ysize = this->zsize = lengthxyz; //these are in world.h

	cells.data.size = lengthxyz / cellsxyz; //.000008m

	/* from here this is a little bit unclear.*/
	grains.data.xPer = 1;//Number of grains in x direction. This is not needed.
	grains.data.yPer = 1;
	grains.data.zPer = 1; // if it were  1 m

	box.xPer = box.yPer = box.zPer = cellsxyz; //these are in grainCA.h
	box.Volume = Vtotal / ncells; //volume per cell (m^3/cell)
	box.size = pow((double)box.Volume, (double)1 / 3.0f);

	long totalcells = box.xPer * box.yPer * box.zPer;

	message << "NumberOfMsGrains:" << grains.data.xPer << "x" << grains.data.yPer << "x" << grains.data.zPer << endl;
	message << "CellsInTotal:" << totalcells << endl;
	message << "CellGridDimensionsXYZ:" << box.xPer << "x" << box.yPer << "x" << box.zPer << endl;
	message << "The automat has a size f:" << this->xsize << "x" << this->ysize << "x" << this->zsize << " m" << endl;
	write(message.str(), &logFile);

	//#as *****added for output form here 

	//ofstream check;
	//check.open("determineCAsizeS.txt");//:as
	//check << "cellsxyz:\t" << cellsxyz << endl;//:as
	//check << "value of Ncells is" << Ncells << endl;//:as
	//check << "ncells:\t" << ncells << endl;//:as
	//check << "Vtotal:\t" << Vtotal << endl;//:as
	//check << "lengthxyz:\t" << lengthxyz << endl;//:as
	//check << "cells.data.size:\t" << cells.data.size << endl;//:as
	//check << "box.xPer:\t" << box.xPer << endl;//:as
	//check << "box.Volume:\t" << box.Volume << endl;//:as
	//check << "box.size:\t" << box.size << endl;//:as
	//check << "total cells:\t" << totalcells << endl;//:as
	//check << "NumberOfMsGrains:" << grains.data.xPer << "x" << grains.data.yPer << "x" << grains.data.zPer << endl;//:as
	//check << "CellsInTotal:" << totalcells << endl;//:as
	//check << "CellGridDimensionsXYZ:" << box.xPer << "x" << box.yPer << "x" << box.zPer << endl;//:as
	//check << "The automat has a size f:" << this->xsize << "x" << this->ysize << "x" << this->zsize << " m" << endl;//:as
	//check.close();//*****added for output till here


}

void caHdl::initializeGrainsSolCA(void) // not clear
{
	grainP ng = new grain(&(this->grains), -1, -1, -1, 0.0, 1e14, 1e14, 1e14, 0, 0, 0, 1, 0.0, 0); // assigned euler angles -1, -1, -1 to powder.
	// (grainPoolP own, Real p1, Real P, Real p2, double rhoGND, double rm0, double ri0, double rw0, int n1, int n2, int n3, short Oriindex, Real NGLS, long DeformGrainIndex )

	grains.defGrainsCount++; // rg: Is this really needed.
	ng->rhoTotal = 1e14;
	ng->CellCount = box.xPer * box.yPer * box.zPer; // new grain cell count is entire domain right now. As we are initalizing the grain
	ng->initialCellNumber = ng->CellCount;
}

void caHdl::startParallelProcesses(int argc, char* argv[])
{
	int mpiErr = MPI_Init(&argc, &argv);

	if (mpiErr != MPI_SUCCESS)
		err->reportError("MPI cannot be initialized");

	int mynode;
	ostringstream message;

	std::time_t time1 = std::time(nullptr);// #added rg
	message << "Starting Parallel Processes " << std::asctime(std::localtime(&time1));// #added rg

	MPI_Comm_size(MPI_COMM_WORLD, &totalNodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

	err->myRank = mynode;

	if (totalNodes == ONEPROCESS)
	{
		if ((boxPer & nodePer & protocol) == (NONBOXPER & NODE_PER & NOCOMM))
		{
			boxPer = BOXPER;
			nodePer = NON_NODE_PER;
			protocol = DETERMINISTIC;
		}
		if (COMM_ENABLED & protocol)
		{
			err->reportWarning("No communication protocol required with only one processor; changed to noComm");
			protocol = NOCOMM;
			transfer.never = NEVER;
		}
		if (boxPer == BOXPER && nodePer == NODE_PER)
		{
			nodePer = NON_NODE_PER;
			err->reportWarning("Only one processor; node periodicity set to box periodicity");
		}
	}

#ifndef CORe_MPI
	if (totalNodes > 1)
		err->reportWarning("Running on more than one processors but no MPI library was compiled. If MPI required, set MPI as a preprocessor directive");
#endif

	if (mynode == MASTER)
	{
		switch (boxPer & nodePer & protocol)
		{
		case BOXPER& NON_NODE_PER& DETERMINISTIC:
			message << "PeriodicityAndProtocol:Periodic box, non periodic nodes and deterministic communication intra-nodes" << endl;
			break;
		case BOXPER& NON_NODE_PER& STATISTICAL:
			message << "PeriodicityAndProtocol:Periodic box, non periodic nodes and statistical communication intra-nodes" << endl;
			break;
		case BOXPER& NON_NODE_PER& NOCOMM:
			message << "PeriodicityAndProtocol:Periodic box, non periodic nodes and no communication intra-nodes: available only for serial simulations" << endl;
			if (totalNodes > 1) err->reportError(ERRTXT("This condition is not possible for parallel simulations"));
			break;
		case NONBOXPER& NODE_PER& NOCOMM:
			message << "PeriodicityAndProtocol:Non-periodic box, periodic nodes and no communication intra-nodes" << endl;
			break;
		case NONBOXPER& NON_NODE_PER& NOCOMM:
			message << "PeriodicityAndProtocol:Non-periodic box, non-periodic nodes and no communication intra-nodes" << endl;
			break;
		case NONBOXPER& NON_NODE_PER& STATISTICAL:
			message << "PeriodicityAndProtocol:Non-periodic box, periodic nodes and statistical communication intra-nodes" << endl;
			break;
		case NONBOXPER& NON_NODE_PER& DETERMINISTIC:
			message << "PeriodicityAndProtocol:Non-periodic box, periodic nodes and deterministic communication intra-nodes" << endl;
			break;
		default:
			write(message.str(), &logFile);
			err->reportError(ERRTXT("No known boundary conditions or communication protocol"));
		}
		write(message.str(), &logFile);
	}

	createMPItypes();
}

void caHdl::createMPItypes(void)
{
	// 1. MPI_IO Oriindex type
	// || MPI_INT			|| MPI_FLOAT				|| MPI_FLOAT			|| MPI_FLOAT
	// || Orientation ID	|| Euler 1					||Euler 2				||Euler 3
	// See io.h struct MPI_IO_OriIndex for more information.
	int elementCounts1[2] = { 1,3 };
	MPI_Aint displacements1[2] = { 0, 4 };
	MPI_Datatype oldTypes1[2] = { MPI_INT, MPI_FLOAT };
	MPI_Type_create_struct(2, elementCounts1, displacements1, oldTypes1, &MPI_IO_OriIndex_Type);


	// 2. MPI_IO Node type // node means process here
	// ||int		||int		||int		||int		||int		||int		||int		
	// ||MPI Rank	||xo		||x_max		||yo		||y_max		||zo		||z_max		
	// See io.h struct MPI_IO_Node for more information.
	int elementCounts2[1] = { 7 };
	MPI_Aint displacements2[1] = { 0 };
	MPI_Datatype oldTypes2[1] = { MPI_INT };
	MPI_Type_create_struct(1, elementCounts2, displacements2, oldTypes2, &MPI_IO_Node_Type);


	// 3. MPI_IO Particle type 
	// ||MPI_LONG		||MPI_LONG		||MPI_LONG		||MPI_LONG		||MPI_LONG		||MPI_LONG		||MPI_LONG		||MPI_LONG
	// it is not utilized
	int elementCountsPar[1] = { 8 };
	MPI_Aint displacementsPar[1] = { 0 };
	MPI_Datatype oldTypesPar[1] = { MPI_LONG };
	MPI_Type_create_struct(1, elementCountsPar, displacementsPar, oldTypesPar, &MPI_IO_Particle_Type);


	// 4. MPI_IO Cell type
	// ||MPI_INT			||MPI_CHAR
	// ||OriIndex			||State(deformend, infected or particle.)
	// See io.h struct MPI_IO_Cell for more info.
	int elementCounts3[2] = { 1,1 };
	MPI_Aint displacements3[2] = { 0, 4 };
	MPI_Datatype oldTypes3[2] = { MPI_INT, MPI_CHAR };
	MPI_Type_create_struct(2, elementCounts3, displacements3, oldTypes3, &MPI_IO_Cell_Type);


	// 5. MPI_IO_Cell Infection state is contiguous
	// ||MPI_FLOAT			||MPI_CHAR
	// ||InfProgress		|| state
	// See io.h MPI_IO_CellInfState for more info.
	int elementCounts4[2] = { 1,1 };
	MPI_Aint displacements4[2] = { 0, 4 };
	MPI_Datatype oldTypes4[2] = { MPI_FLOAT, MPI_CHAR };
	MPI_Type_create_struct(2, elementCounts4, displacements4, oldTypes4, &MPI_IO_CellInfState_Type);


	// 6. MPI_IO_Grain type //###LB: New MPI_type for the transfer of grains
	// see io.h struct MPI_IO_GrainInfo for more info.
	int elementCountsGrain[4] = { 17,7,5,7 };

	MPI_Aint displacementsGrain[4] = { offsetof(MPI_IO_GrainInfo, phi1_TB), offsetof(MPI_IO_GrainInfo, Oriindex),
									   offsetof(MPI_IO_GrainInfo, alreadyNuc), offsetof(MPI_IO_GrainInfo, oriIndex) };

	MPI_Datatype oldTypesGrain[4] = { MPI_DOUBLE, MPI_SHORT, MPI_LONG, MPI_INT };
	MPI_Type_create_struct(4, elementCountsGrain, displacementsGrain, oldTypesGrain, &MPI_IO_GrainInfo_Type);


	// 7. MPI_IO_Nucleus
	// || MPI_SHORT		|| MPI_SHORT		|| MPI_SHORT		||MPI_LONG
	// ||x				||y					||z					||ig(dont know what is this)
	// See grainCA.h struct exportNucleus for more info.
	int elementCountsNuc[2] = { 3, 1 };
	MPI_Aint displacementNuc[2] = { offsetof(exportNucleus, x),offsetof(exportNucleus, ig) };
	MPI_Datatype oldTypesNuc[2] = { MPI_SHORT, MPI_LONG };
	MPI_Type_create_struct(2, elementCountsNuc, displacementNuc, oldTypesNuc, &MPI_IO_Nucleus_Type);

	//8. MPI_IO_Ori
	// ||MPI_DOUBLE		||MPI_DOUBLE		||MPI_DOUBLE		||MPI_DOUBLE		||MPI_DOUBLE		||MPI_DOUBLE		||MPI_DOUBLE		||MPI_LONG		
	// ||phi1			|| PHI				||phi1				||q0()quaternian	||q1				||q2				||q3				||count	
	// See grainCA.h struct exportOri for more info.
	int elementCountsOri[2] = { 7, 1 };
	MPI_Aint displacementOri[2] = { offsetof(exportOri, phi1),offsetof(exportOri, count) };
	MPI_Datatype oldTypesOri[2] = { MPI_DOUBLE, MPI_LONG };
	MPI_Type_create_struct(2, elementCountsOri, displacementOri, oldTypesOri, &MPI_IO_Ori_Type);

	// Now here we commit all MPI data types.

	// writing to mpifile data types
	MPI_Type_commit(&MPI_IO_OriIndex_Type); // for writing data
	MPI_Type_commit(&MPI_IO_Node_Type); // for writing data
	MPI_Type_commit(&MPI_IO_Particle_Type); // not used.
	MPI_Type_commit(&MPI_IO_Cell_Type);// for writing data
	MPI_Type_commit(&MPI_IO_CellInfState_Type); //  for writing data, Cell infection state.


	//reading from mpifile datatypes.
	MPI_Type_commit(&MPI_IO_GrainInfo_Type);// for reading data from binary file
	MPI_Type_commit(&MPI_IO_Nucleus_Type);// for reading data from binary file
	MPI_Type_commit(&MPI_IO_Ori_Type); // for reading data from binary file
}

void caHdl::prepareAutomata(void)
{
	ostringstream message;

	message << "Determing processors lattice" << endl;
	determineNodeTopology(); // not clear

	message << "Initial Temperature" << endl;
	initializeTemperatureCells();

	message << "Initializing cells" << endl;
	initializeCells();

	message << "Inserting Nuclei" << endl;
	insertNucleiSolidification();

	message << "Initial Texture" << endl;
	determineInitialTexture();
	message << "Porosity checking" << endl; //as: added
	//InsufficientMelting();
	
	if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));

	write(message.str(), &logFile);
}

void caHdl::determineNodeTopology(void) //?? not clear
{
	//******************
	/*ofstream check;
	check.open("determineNodeTopology.txt");*/
	//***********************
	ostringstream message;

	int np = totalNodes; // total number of processes  
	int mynode; // current process rank.

	//*****
	//check << "total nodes are:" << np <<"\n\n"<< endl;
	//check << "my node is:" << mynode << "\n\n" << endl;
	//*******



	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

	long partition[3] = { 0 };
	factorizeIn3((long)totalNodes, partition);

	parallelHdl = new processHdl(this, np, partition[0], partition[1], partition[2], box.xPer, box.yPer, box.zPer);

	int nodeMesh = parallelHdl->xPer * parallelHdl->yPer * parallelHdl->zPer; // Number of total processes 1

	//***** #as added tp read 
	//check << "nodeMesh are:" << nodeMesh << "\n\n" << endl; //***** #as added to read 

	//*******

	parallelHdl->all = (parallelDataP*)malloc(nodeMesh * sizeof(parallelDataP));

	if (!parallelHdl->all) err->reportError(ERRTXT("Cannot allocate more memory"));

	Real xtol = box.xPer * TOLFP(box.xPer); //as:499 why ?? 
	Real ytol = box.yPer * TOLFP(box.yPer);
	Real ztol = box.zPer * TOLFP(box.zPer);

	Real xinc = ((Real)box.xPer / partition[0]); //as:number of cells per subdomain along the x-edge 500
	Real yinc = ((Real)box.yPer / partition[1]);
	Real zinc = ((Real)box.zPer / partition[2]);
	//*****
	/*check << "TOLFP is:\t" << TOLFP(box.xPer) << "\n\n"<<endl;//***** #as added to read
	check << "nodeMesh are:" << nodeMesh << "\n\n" << endl;//***** #as added to read
	check << "xtol  ytol  ztol are::\t" << xtol << "\t" << ytol << "\t" << ztol << "\n\n"<< endl; //***** #as added to read
	check << "xinc  yinc and zinc are::\t" << xinc << "\t" << yinc << "\t" << zinc<<"\n\n" << endl; //***** #as added to read
	check << "partiions values are\t" << partition[0] << "\t" << partition[1] << "\t" << partition[2] << "\n" << endl;*/ //***** #as added to read 
	//*****


	int k = 0;
	for (Real z = 0; z < ztol; z += zinc) {
		for (Real y = 0; y < ytol; y += yinc) {
			for (Real x = 0; x < xtol; x += xinc)
			{// we use this loop to define from the cells within each processes.
				long minx = (long)x;
				long maxx = (long)(x + xinc);
				long miny = (long)y;
				long maxy = (long)(y + yinc);
				long minz = (long)z;
				long maxz = (long)(z + zinc);

				int nix = (int)(x / xinc); //equals a node number in x-direction 0 <= nix (<= nx)? there are max. 255 nodes in one direction
				int niy = (int)(y / yinc);
				int niz = (int)(z / zinc);

				parallelHdl->all[k] = new parallelData(minx, maxx, miny, maxy, minz, maxz, k, nix, niy, niz);

				//Use of memset is only possible because xnodes, ynodes and znodes are of type char; max. 255x255x255 nodes. Check for portability ??
				if (y == 0 && z == 0)	memset(&(parallelHdl->xnodes[minx]), (unsigned char)(x / xinc), (maxx - minx) * sizeof(unsigned char));
				if (x == 0 && z == 0)	memset(&(parallelHdl->ynodes[miny]), (unsigned char)(y / yinc), (maxy - miny) * sizeof(unsigned char));
				if (x == 0 && y == 0)	memset(&(parallelHdl->znodes[minz]), (unsigned char)(z / zinc), (maxz - minz) * sizeof(unsigned char));

				//six neighboring nodes towards the <100> subdomain faces
				parallelHdl->all[k]->neighbourNodes[C_RIGHT] = parallelHdl->getNextNode(nix, niy, niz, 1, 0, 0); //stores the node-id
				parallelHdl->all[k]->neighbourNodes[C_LEFT] = parallelHdl->getNextNode(nix, niy, niz, -1, 0, 0);
				parallelHdl->all[k]->neighbourNodes[C_TOP] = parallelHdl->getNextNode(nix, niy, niz, 0, 1, 0);
				parallelHdl->all[k]->neighbourNodes[C_BOTTOM] = parallelHdl->getNextNode(nix, niy, niz, 0, -1, 0);
				parallelHdl->all[k]->neighbourNodes[C_FRONT] = parallelHdl->getNextNode(nix, niy, niz, 0, 0, -1);
				parallelHdl->all[k]->neighbourNodes[C_REAR] = parallelHdl->getNextNode(nix, niy, niz, 0, 0, 1);

				if (mynode == k) {

					this->x0 = minx; // these values are in world.h
					this->y0 = miny;
					this->z0 = minz;
					this->xPer = maxx - minx;
					this->yPer = maxy - miny;
					this->zPer = maxz - minz;
					this->zFac = this->xPer * this->yPer;
					this->rank = mynode;

					// Now if we have more than one processes and these proceses communicate data with each other throuch interface.
					if (totalNodes > ONEPROCESS && (protocol & COMM_ENABLED))
					{
						//six neighboring nodes towards the <100> subdomain faces, 2D oriindex array
						nodeFace[C_TOP] = new worldInterface(this, this->xPer, this->zPer, C_TOP);
						nodeFace[C_BOTTOM] = new worldInterface(this, this->xPer, this->zPer, C_BOTTOM);

						nodeFace[C_LEFT] = new worldInterface(this, this->yPer, this->zPer, C_LEFT);
						nodeFace[C_RIGHT] = new worldInterface(this, this->yPer, this->zPer, C_RIGHT);

						nodeFace[C_FRONT] = new worldInterface(this, this->xPer, this->yPer, C_FRONT);
						nodeFace[C_REAR] = new worldInterface(this, this->xPer, this->yPer, C_REAR);

						nodeFace[C_RIGHT]->nextNode = parallelHdl->all[k]->neighbourNodes[C_RIGHT];
						nodeFace[C_LEFT]->nextNode = parallelHdl->all[k]->neighbourNodes[C_LEFT];

						nodeFace[C_TOP]->nextNode = parallelHdl->all[k]->neighbourNodes[C_TOP];
						nodeFace[C_BOTTOM]->nextNode = parallelHdl->all[k]->neighbourNodes[C_BOTTOM];

						nodeFace[C_FRONT]->nextNode = parallelHdl->all[k]->neighbourNodes[C_FRONT];
						nodeFace[C_REAR]->nextNode = parallelHdl->all[k]->neighbourNodes[C_REAR];
					}
				}
				k++;
			}
		}
	}

	if (mynode == MASTER)
	{
		if (partition[0] > 4 && partition[1] == 1 && partition[2] == 1) message << "Prime number of processors is not efficient\n";
		message << "Automaton subdivided in a:" << partition[0] << "x" << partition[1] << "x" << partition[2] << " processors lattice\n";
		write(message.str(), &logFile);
	}
}

void caHdl::initializeTemperatureCells(void)
{
	////******************** as added on 2/12/2022
	//	ofstream check;// as added for ouput
	//	check.open("initializeTemperatureCells.txt");// created a file with this name
	//	if (!check) {  // to check it worked/
	//		cout << "File not created!";// as added for ouput
	//	}
	//	else { cout << "successsinitializeTemperatureCells "; }// as added for ouput
	// //added on 12/1/22
	////************************added on 2/12/2022


	ostringstream message;

	message << "Initializing temperature cells" << endl;

	long ncells = this->Tcells->ntcells;// number of temperature cells  3000 20*15*10
	long nnodes = this->Tcells->nodes->size;// number of temperature nodes 3696 21*16*11

	//************************added on 2/12/2022
	//check << "ncells are" << setw(12) << ncells << setw(12) << "nnodes are:" << setw(12) << nnodes << "\n" << endl; //added on 12/1/22

	//************************added on 2/12/2022

	Real xsize = Tcells->xmax - Tcells->xmin; // total domain size in actual frames
	Real ysize = Tcells->ymax - Tcells->ymin;
	Real zsize = Tcells->zmax - Tcells->zmin;

	long xcells = box.xPer;// number of automaton cells in x direction 500
	long ycells = box.yPer;
	long zcells = box.zPer;

	Real factorx = xcells / xsize;// how many automaton cells if we have unit length of automaton domain 1250000
	Real factory = ycells / ysize;
	Real factorz = zcells / zsize;

	//**********//added on 12/1/22 as

	/*	check << "x size y size and z size are" << xsize << setw(10) << ysize << setw(10) << zsize << "\n" << endl;//added on 12/1/22
		check << "xcells  ycelss and   zcells are:"<< xcells << setw(10) << ycells << setw(10) << zcells << setw(10) << "\n" << endl;//added on 12/1/22
		check << "factors are:" << setw(10) << factorx << setw(10) <<"\n"<< endl;
		check << "\n  value of x          value of y          valuye of z        xi        yi       zi        \n" << endl;
	*/
	//************** as	//added on 12/1/22


	this->Tcells->all = (tempCellP*)calloc(ncells, sizeof(tempCellP));// assigns memory for temperature cell pointers
	this->Tcells->nodes->all = (TNodeP*)calloc(nnodes, sizeof(TNodeP));// assigns memory for temperature node pointers

	if (!(this->Tcells->all)) err->reportError(ERRTXT("Cannot allocate more memory"));
	if (!(this->Tcells->nodes->all)) err->reportError(ERRTXT("Cannot allocate more memory"));

	for (int i = 0; i < nnodes; i++)// for all the temperature nodes.
	{
		TNodeRawData rawData = this->Tcells->nodes->tnodesraw[i];

		Real x = rawData.pos.x - Tcells->xmin;
		Real y = rawData.pos.y - Tcells->ymin; //ymin is .-0004 (y-ymin)
		Real z = rawData.pos.z - Tcells->zmin;

		long xi = (long)((x * factorx) + 0.5);//factorx 1250000
		long yi = (long)((y * factory) + 0.5);
		long zi = (long)((z * factorz) + 0.5);

		x = (Real)xi;
		y = (Real)yi;
		z = (Real)zi;

		if (xi > 0 && xi < xcells) x = ((Real)(xi - 1)) + 0.5;
		if (yi > 0 && yi < ycells) y = ((Real)(yi - 1)) + 0.5;
		if (zi > 0 && zi < zcells) z = ((Real)(zi - 1)) + 0.5;

		QUICKASSERT(i == rawData.id);

		TNodeP Tnode = new TNode(x, y, z, rawData.id);
		this->Tcells->nodes->all[i] = Tnode;

		//*********** added on 1/12/2022
		//check << x << setw(10) << y << setw(10) << z << setw(10) <<  xi << setw(10) << yi << setw(10) << zi << setw(10) << "\n"<< endl;
		//*********** added on 1/12/2022


	}



	short* xdi = (short*)calloc((size_t)(Tcells->xPerT + (size_t)1), sizeof(short));
	short* ydi = (short*)calloc((size_t)(Tcells->yPerT + (size_t)1), sizeof(short));
	short* zdi = (short*)calloc((size_t)(Tcells->zPerT + (size_t)1), sizeof(short));

	// allocates which temperature node lies in which automaton cell. I guess.
	for (int i = 0; i <= Tcells->xPerT; i++) {
		xdi[i] = (short)(((Tcells->xdiscr[i] - Tcells->xmin) * factorx) + 0.5);

	}


	for (int i = 0; i <= Tcells->yPerT; i++)
	{
		ydi[i] = (short)(((Tcells->ydiscr[i] - Tcells->ymin) * factory) + 0.5);



	}

	for (int i = 0; i <= Tcells->zPerT; i++)
	{
		zdi[i] = (short)(((Tcells->zdiscr[i] - Tcells->zmin) * factorz) + 0.5);

	}


	short counter = 0;

	//**********added on 1/12/2022
	ofstream check3;
	check3.open("temperaturecells check3.txt");
	if (!check3)
		cout << "not created\n" << endl;

	check3 << "Tcells->xPerT and Tcells->yPerT and Tcells->zPerT are :\n" << Tcells->xPerT << setw(15) << Tcells->yPerT << setw(15) << Tcells->zPerT << setw(15) << "\n" << endl;
	//**********added on 1/12/2022

	for (int zi = 0; zi < Tcells->zPerT; zi++) // for all the temperature cells
	{
		//check3 << "for zi= (" << zi << ")\n" << endl;// as:ouput check

		for (int yi = 0; yi < Tcells->yPerT; yi++)
		{
			//check3 << " and for yi= (" << yi << ")\n" << endl;// as:ouput check

			for (int xi = 0; xi < Tcells->xPerT; xi++)
			{
				//check3 << " and for xi= (" << xi << ")\n" << endl; // as:ouput check

				tempCellP Tcell = new tempCell(Tcells, xdi[xi], xdi[xi + 1] - 1, ydi[yi], ydi[yi + 1] - 1, zdi[zi], zdi[zi + 1] - 1);

				cellNodesRawData rawCellData = Tcells->tcellraw[counter]; // neighbouring temperature nodes of the given temperature cell 

				Tcell->id = counter;

				short parallelNodes[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };

				char maskx = 0x4;
				char masky = 0x2;
				char maskz = 0x1;

				for (short i = c000; i <= c111; i++)// for all the neighbouring nodes of the current temperature cell c000 is 1 and so on (not clear)
				{
					char displx = 0x0;
					char disply = 0x0;
					char displz = 0x0;

					displx = (i & maskx) / maskx; // ??
					disply = (i & masky) / masky;
					displz = (i & maskz) / maskz;

					long corner = rawCellData.cellnodes[i];
					Tcell->data.vertexTemp[i] = Tcells->nodes->all[corner];

					short xnode = xdi[xi + displx] - displx;
					short ynode = ydi[yi + disply] - disply;
					short znode = zdi[zi + displz] - displz;

					int noderank = parallelHdl->cellNode((long)xnode, (long)ynode, (long)znode); // find process in which the node lies

					parallelNodes[i] = noderank;
					//check3 << "parallel nodes(" << i << ") is:" << noderank << "\n";//added for output check
					//check3 << "xnode,y,z" << xnode << " " << ynode << "  " << znode <<"\n"<< endl; //added for output check
				}

				if (isMyRankInTheList(parallelNodes, 8))
				{ // if the temperature node lies in that process, give the data corresponding of that temp-node to the process.
					// add that temperature cell the bucket.
					Tcell->allocateTemperatureTimeLocalNodes(Tcells->T_timer.size);
					Tcells->addTempCelltoBucket(Tcell);
				}

				QUICKASSERT((short)ceil(Tcell->data.vertexTemp[c000]->pos.x) == xdi[xi]);
				QUICKASSERT((short)floor(Tcell->data.vertexTemp[c100]->pos.x) == xdi[xi + 1] - 1 || xdi[xi + 1] == xcells);
				QUICKASSERT((short)ceil(Tcell->data.vertexTemp[c000]->pos.y) == ydi[yi]);
				QUICKASSERT((short)floor(Tcell->data.vertexTemp[c010]->pos.y) == ydi[yi + 1] - 1 || ydi[yi + 1] == ycells);
				QUICKASSERT((short)ceil(Tcell->data.vertexTemp[c000]->pos.z) == zdi[zi]);
				QUICKASSERT((short)floor(Tcell->data.vertexTemp[c001]->pos.z) == zdi[zi + 1] - 1 || zdi[zi + 1] == zcells);

				Tcells->all[counter] = Tcell;

				counter++;
			}
		}
	}

	free(Tcells->tcellraw);
	Tcells->tcellraw = NULL;
	free(Tcells->nodes->tnodesraw);
	Tcells->nodes->tnodesraw = NULL;
	free(xdi);
	free(ydi);
	free(zdi);
	free(this->Tcells->nodes->temperatureSeriesRawData);
	write(message.str(), &logFile);
}



bool caHdl::isMyRankInTheList(short* list, int size) //not clear
{
	//ofstream check4; //  as: to check the output
	//check4.open("isMyRankInTheList.txt");// as: to check the output
	int myRank = myNode();
	for (int i = 0; i < size; i++)
	{
		//check4 << "my rank value is" << list[i] << "\n" << endl;// as: to check the output
		if (myRank == list[i])
			return true;



	}
	return false;

}

void caHdl::initializeCells(void)
{
	//// added to see output not required
/*	ofstream check5; //as: added to see output not required
	check5.open("initializeCellscheck5.txt");//as: added to see output not required
	*///if(!check5) cout << "not created\n" << endl;//as: added to see output not required


	ostringstream message;
	message << "CA cells initialization" << endl;

	long partitionCells = 0;
	long i;

	// check5 << "TIMESLOTCOUNT is" << TIMESLOTCOUNT << "\n" << endl;// added to se eoutput

	for (i = 0; i <= TIMESLOTCOUNT; i++) { cells.first[i] = NULL; }
	cells.data.count = 0;  // cells is struct of AutomatonCellPool, mind it, it is not AutomatonCell.
	cells.data.retiredCount = 0;
	cells.firstRecycled = NULL;
	cells.countRecycled = 0;

	cells.data.xPer = this->xPer;//gives the number of CA cells in my subdomain/ process.

	//check5 << "cells.data.xper" << cells.data.xPer << "\n" << endl; //added to see output

	cells.data.yPer = this->yPer;

	//check5 << "cells.data.yper" << cells.data.yPer << "\n" << endl;//added to see output

	cells.data.zPer = this->zPer;

	//check5 << "cells.data.zper" << cells.data.zPer << "\n" << endl;//added to see output

	partitionCells = cells.data.xPer * cells.data.yPer * cells.data.zPer; // total number of CA cells in the entire subdomain which is 12500000.

	//check5 << "partitionCells are" << partitionCells << "\n"; // added to see output

	cells.all = (cellP*)malloc(partitionCells * sizeof(cellP)); //cells as a VRM continuous chunk of 4Byte cellP pointers

	if (!cells.all) err->reportError(ERRTXT("Cannot allocate more memory"));

	cells.data.zFac = cells.data.xPer * cells.data.yPer;

	grains.data.xCellCount = box.xPer; // this should not be used. Was basically meant for crystallization simulation.
	grains.data.yCellCount = box.yPer;
	grains.data.zCellCount = box.zPer;

	// check5 << "Box.xPer are: " << box.xPer << endl;  //added to check the output

	message << "CA cells temperature assignation" << endl;
	establishCellOwnershipToTemperatureCell();

	write(message.str(), &logFile);
}

void caHdl::establishCellOwnershipToTemperatureCell(void)
{	
	long count = 0;
	for (int i = 0; i < Tcells->PNBucket.size; i++) // for all the temperature cells in your process/subdomain 
	{

		tempCellP tc = Tcells->PNBucket.TCparallelNode[i]; //

		short xl = tc->mySpace.xleft; // these are the beginning and the end of the temperature_cells in all direction. not clear
		short xr = tc->mySpace.xright;
		short yl = tc->mySpace.yleft;
		short yr = tc->mySpace.yright;
		short zl = tc->mySpace.zleft;
		short zr = tc->mySpace.zright;

		
		uint32_t myTCell = (uint32_t)tc->id;//id was short i.e 16 bits here it is convereted in 32 bits.

		short ixavg = (xl + (long)(0.5 * (xr - xl)));// finding the midpoint of the temperature cell.
		short iyavg = (yl + (long)(0.5 * (yr - yl)));
		short izavg = (zl + (long)(0.5 * (zr - zl)));


		Real avgTimeToMelt = tc->determineNucleationTime(ixavg, iyavg, izavg); // refer tempControl.cpp

		
		uint_fast32_t timeDiscrete = 0x3FFF; // 16383

		if ((long)avgTimeToMelt != LIQUID && (long)avgTimeToMelt != NO_MELT )
			timeDiscrete = (uint_fast32_t)(avgTimeToMelt * Tcells->timeEncodeFactor);

		timeDiscrete <<= 1; // Multiplying timeDiscrete by 2 or left shift by 1 bit.

		QUICKASSERT(Tcells->all[myTCell] == tc);

		//*********** //added for output check on 7/12/2022
		//ofstream check;
		//	check.open("establishCellOwnershipToTemperatureCellpart2.txt");
		//	check << "x\t\t\t y\t\t\t z\t\t\t index\t\t\t testTemp " << endl; //added for output check on 7/12/2022
			////****************

		for (short iz = zl; iz <= zr; iz++)// for all the ca cell in my process
		{


			for (short iy = yl; iy <= yr; iy++)
			{

				for (short ix = xl; ix <= xr; ix++)
				{

					long index = (iz - z0) * cells.data.zFac + (iy - y0) * cells.data.xPer + (ix - x0);// find the local index of the ca cell.

					//******** //added for output check on 7/12/2022
					// check << ix << "\t\t\t" << iy << "\t\t\t" << iz << "\t\t\t" << index<<"\t\t\t";//added for output check on 7/12/2022
					//*************//added for output check on 7/12/2022					
					int cellnode = parallelHdl->cellNode(ix, iy, iz);// find in which process that ca cell lies.

					if (cellnode != myNode()) continue;

					cells.all[index] = (cellP)0x00000001;

					uint32_t TCellIdxShift = myTCell << 15;

					uint32_t ci = (uint32_t)cells.all[index];

					QUICKASSERT(isCellSolid(cells.all[index]));

					ci |= TCellIdxShift;
					ci |= timeDiscrete;

					// ||xxxx|xxxx|xxxx|xxxx|x yyy|yyyy|yyyy|yyy|1||
					//  Temperature Cell Id   |  timeDiscrete     |last bit(==1) showing initialized ca Cell belongs to deformeg grain.
					cells.all[index] = (cellP)ci;

					Real testTemp = Tcells->getCellNucleationTime((cellP)ci);//testTemp never used.
//					this->PT.tmax.push_back(testTemp);
					//cout << "testTemp:" << testTemp << endl;
					//cout << "test2" << t << endl;
					//check << testTemp << endl;
					count++;
					//cout << cells.all[index]->T << endl;
				}

			}
		}
	}

	QUICKASSERT(count == xPer * yPer * zPer);
}

void caHdl::insertNucleiSolidification(void)
{
	
	//InsufficientMelting();
	insertBulkNucleiSolidification();
	//insertSurfaceNucleiSolidification();
	//insertEpitaxialNucleationSolidification();
}


//probablisitic nucleation layer by layer with grad

//void caHdl::insertBulkNucleiSolidification(void)    
//{
//#define NUMBER_HORIZONTAL_PASSES  5
//#define NUMBER_VERTICLE_PASSES  12//12 //for single layer it is 1
//#define SPACING_BETWEEN_HORIZONTAL_PASSES 70/0.8
//#define SPACING_BETWEEN_VERTICLE_PASSES 30/0.8
//#define LASER_START_X 50/0.8
//#define LASER_START_Y 70/0.8
//	ofstream nucleusData("NucleusData.txt", ios::out);
//	long bulkNucNumber = 23000/2;// 23000// number of nuclei
//	Real shortestNucTime = 1e10;//??
//	Real prob=0;
//	Real ptest=0;
//	long minNucx = 1000;
//	long maxNucx = 0;  //WHY max is 0 ??
//	long minNucy = 1000;
//	long maxNucy = 0;
//	long minNucz = 1000;
//	long maxNucz = 0;
//	int n = 1;
//	Real tmax = 0;
//	Real maxtime = Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]; //as:for pob nuc
//	Real tt = 0;// maxtime;//as for prob nuc
//	//Real cX, rx, ry, rz, dY, eZ, fz, phia, phib, phic, meltime; //**as added for output check
//	//cellP cz = NULL; //** as added for output check
//	//TCellP T=NULL;  //asaded for prb nuc
//	long key;//as:added for prb nuc
//	long surfaceNucNumber = 23000/2 ;// total number
//	long bulkNucleation = 23000/2 ;
//	long numberNucleusEachLayer = bulkNucleation / NUMBER_VERTICLE_PASSES;
//	long numberofNucleusEachLine = numberNucleusEachLayer / NUMBER_HORIZONTAL_PASSES;
//	vector <int> y;
//	vector <int> x;
//
//
//
//	for (int i = 0; i < NUMBER_VERTICLE_PASSES; i++) {
//		y.push_back(LASER_START_Y + i * SPACING_BETWEEN_VERTICLE_PASSES);
//		for (int j = 1; j < NUMBER_HORIZONTAL_PASSES+1; j++) {
//			x.push_back(LASER_START_X + j * SPACING_BETWEEN_HORIZONTAL_PASSES);
//
//		}
//
//	}
//
//
//	while (bulkNucNumber > 0) // ( TimeT < maxtime) //(bulkNucNumber > 0)
//	{
//		short success = 0;
//		do {
//			if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));
//			int indexY = int((bulkNucleation - bulkNucNumber) / numberNucleusEachLayer); // this decides the maximum numbe of nucleus per layer
//			if (indexY < 0) indexY = 0;
//			if (indexY >= NUMBER_VERTICLE_PASSES) indexY = NUMBER_VERTICLE_PASSES - 1;
//			int indexX = int((bulkNucleation - bulkNucNumber) / numberofNucleusEachLine); // check it works ?this decides the maximum numbe of nucleus per layer
//			if (indexX < 0) indexX = 0;
//			if (indexX >= NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES) indexX = NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES - 1;
//			QUICKASSERT(indexX <= NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES);
//			
//			long ix =  (long)(x[indexX]-this->r.parkMiller() *(SPACING_BETWEEN_HORIZONTAL_PASSES ));//parkmiller random number 
//			
//			
//			
//
//			/*
//			long ix = (long)(this->r.parkMiller() * box.xPer);//parkmiller random number 
//			long iz = (long)(this->r.parkMiller() * box.zPer);
//			long iy = (long)(this->r.parkMiller() * box.yPer);
//			*/
//
//
//			QUICKASSERT(ix >= 0 && ix < box.xPer);
//	//	QUICKASSERT(iz >= 0 && iz < box.zPer);
//
//			//int cellnode = parallelHdl->cellNode(ix, iy, iz); // find in which process that cell lies. WHICH CELL?
//
//			int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());// converting decimal to whole number smaller than the number itself
//
//			
//
//
//			Real phi1_fixed = grainOriNuc[indexOri + 0];
//			Real PHI_fixed = grainOriNuc[indexOri + 1];
//			Real phi2_fixed = grainOriNuc[indexOri + 2];
//
//			for (int i = 0; i < y.size(); i++)
//			{
//				long iz = (long)(this->r.parkMiller() * (box.zPer));
//				long iy = y[indexY]-(long)(this->r.parkMiller() * (SPACING_BETWEEN_VERTICLE_PASSES));
//				if (iy >= box.yPer) {
//					iy = iy - SPACING_BETWEEN_VERTICLE_PASSES;
//				}
//				QUICKASSERT(iy >= 0 && iy < box.yPer);
//				int cellnode = parallelHdl->cellNode(ix, iy, iz); // find in which process that cell lies.
//
//				Real phi1 = phi1_fixed+(float)rand() / RAND_MAX;
//				Real PHI = PHI_fixed+(float)rand() / RAND_MAX;
//				Real phi2 = phi2_fixed +(float)rand() / RAND_MAX;
//
//
//				long DGI = -1; //Not needed for solidification
//
//				int oi = 1; //*Index represents the affiliation to an ideal ori*/
//
//				MPI_IO_GrainInfo grainData;// this is a struct.
//
//				/*ix = 203;
//				iy = 417;
//				iz = 235;*/
//
//
//
//				if (myNode() == cellnode) // if the generated coodinates of the nucleus lies in my process
//				{
//					cellP* cp = getCellPerP(ix, iy, iz);
//					cellP c = *cp;
//					cellP newCell = NULL;
//					if (!cellExists(c))
//					{
//						tempCellP TCell = getTCellofCell(c);
//
//						Real timeAtSolidification = getNucTimeInactiveCell(c);
//
//
//
//						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]) //as: what is it doing ?
//						{
//
//							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
//							// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)
//
//
//							//#as in tetsing: prob.nucleation
//
//							long len = Tcells->T_timer.size;
//
//							for (long i = 0; i < len - 1; i++)
//							{
//								if ((.0000016 + timeAtSolidification) >= Tcells->T_timer.timeSeries[i] && (.0000016 + timeAtSolidification) <= Tcells->T_timer.timeSeries[i + 1])
//								{
//									key = i;
//									break;
//								}
//							}
//							long prevKey = 0;
//							long nextKey = key + 1;
//
//							if (key == Tcells->T_timer.size - 1) nextKey = key;
//
//							Real ct = TCell->interpolateTemperatureSpace(ix, iy, iz, (.0000016 + timeAtSolidification), key, nextKey);  //#as in tetsing: prob.nucleation
//							Real undercooling = gen.Tmelt - ct;
//							//prob = .85;
//							prob = this->nucprobability.calculateCumulativeProbability(undercooling, 5, 1);
//
//							//Real prob1 = (100 * (static_cast<double>(1) / 1 * sqrt(2 * 22 / 7)) * exp(-(SQR((1375 - ct) - 5) / (2 * SQR(1))))); // gaussian
//
//							ptest = (float)rand() / RAND_MAX;       //random number between 1 and 0
//
//							if (prob > ptest)
//							{
//
//								if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;
//
//								if (minNucx > ix) minNucx = ix;
//								if (minNucy > iy) minNucy = iy;
//								if (minNucz > iz) minNucz = iz;
//								if (maxNucx < ix) maxNucx = ix;
//								if (maxNucy < iy) maxNucy = iy;
//								if (maxNucz < iz) maxNucz = iz;
//
//
//
//								if (newCell)
//								{
//									nucleusData << 0.8 * ix << "\t" << 0.8 * iy << "\t" << 0.8 * iz << "\t" << phi1 << "\t" << PHI << "\t" << phi2 << endl; //rg.: Multiplication by 0.8 to convet to abbsolute coordiantes.
//									success = 1;
//									grains.Nuc++;
//									prepareGrainForBroadcast(newCell, &grainData);
//									bulkNucNumber--;
//									
//								}
//
//							}
//
//
//						}
//					}
//				}
//
//				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
//				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);
//
//				if (myNode() != cellnode && success)
//				{
//					addGrain(&grainData);
//					grains.Nuc++;
//				}
//			}
//		} while (!success);
//				
//		//bulkNucNumber--;
//	
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	Real reducedTime = 0;
//
//	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//
//	time = reducedTime;
//	Tcells->setTtimer(shortestNucTime);
//	nucleusData.close(); // #rg
//	write_nucleus_data_SEPI();
//}

//probablisitic nucleation layer by layer without grad

void caHdl::insertBulkNucleiSolidification(void)
{
#define NUMBER_HORIZONTAL_PASSES  5
#define NUMBER_VERTICLE_PASSES  12//12 //for single layer it is 1
#define SPACING_BETWEEN_HORIZONTAL_PASSES 70/0.8
#define SPACING_BETWEEN_VERTICLE_PASSES 30/0.8
#define LASER_START_X 72//57/.8
#define LASER_START_Y 70/0.8
#define Meltpool_X 80/.8
	ofstream nucleusData("NucleusData.txt", ios::out);
	long bulkNucNumber = 2*23000;// 23000// number of nuclei
	Real shortestNucTime = 1e10;//??
	Real prob = 0;
	Real ptest = 0;
	long minNucx = 1000;
	long maxNucx = 0;  //WHY max is 0 ??
	long minNucy = 1000;
	long maxNucy = 0;
	long minNucz = 1000;
	long maxNucz = 0;
	int n = 1;
	Real tmax = 0;
	Real maxtime = Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]; //as:for pob nuc
	Real tt = 0;// maxtime;//as for prob nuc
	//Real cX, rx, ry, rz, dY, eZ, fz, phia, phib, phic, meltime; //**as added for output check
	//cellP cz = NULL; //** as added for output check
	//TCellP T=NULL;  //asaded for prb nuc
	long key;//as:added for prb nuc
	long surfaceNucNumber = 2*23000;// total number
	long bulkNucleation = 2*23000;
	long numberNucleusEachLayer = bulkNucleation / NUMBER_VERTICLE_PASSES;
	long numberofNucleusEachLine = numberNucleusEachLayer / NUMBER_HORIZONTAL_PASSES;
	vector <int> y;
	vector <int> x;



	for (int i = 0; i < NUMBER_VERTICLE_PASSES; i++) {
		y.push_back(LASER_START_Y + i * SPACING_BETWEEN_VERTICLE_PASSES);
		for (int j =0; j < NUMBER_HORIZONTAL_PASSES ; j++) {
			x.push_back(LASER_START_X + j * SPACING_BETWEEN_HORIZONTAL_PASSES);

		}

	}


	while (bulkNucNumber > 0) // ( TimeT < maxtime) //(bulkNucNumber > 0)
	{
		short success = 0;
		do {
			if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));
			int indexY = int((bulkNucleation - bulkNucNumber) / numberNucleusEachLayer); // this decides the maximum numbe of nucleus per layer
			if (indexY < 0) indexY = 0;
			if (indexY >= NUMBER_VERTICLE_PASSES) indexY = NUMBER_VERTICLE_PASSES - 1;
			int indexX = int((bulkNucleation - bulkNucNumber) / numberofNucleusEachLine); // check it works ?this decides the maximum numbe of nucleus per layer
			if (indexX < 0) indexX = 0;
			if (indexX >= NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES) indexX = NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES - 1;
			QUICKASSERT(indexX <= NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES);
			//long ix = 0;
			/*if (x[indexX] == LASER_START_X)
			{
				 ix = (long)(this->r.parkMiller() * Meltpool_X);
			}
			else
			{
				ix = (long)(x[indexX] - .5 * Meltpool_X) + (long)(this->r.parkMiller() * (Meltpool_X));//parkmiller random number


			} */
			 
			long iy = (long)(y[indexY] - ((long)(this->r.parkMiller() * (SPACING_BETWEEN_VERTICLE_PASSES))));
			long ix = (long)(this->r.parkMiller() * box.xPer);
			if (iy >= box.yPer)
				{
				iy = iy - SPACING_BETWEEN_VERTICLE_PASSES;
				}
	
				if (iy <0 )
				{
				 iy=0;
				}


			long iz = (long)(this->r.parkMiller() * (box.zPer));
			
			QUICKASSERT(iy >= 0 && iy < box.xPer);
			QUICKASSERT(ix >= 0 && ix < box.xPer);
			QUICKASSERT(iz >= 0 && iz < box.zPer);
			
			
			int cellnode = parallelHdl->cellNode(ix, iy, iz); // find in which process that cell lies.
			int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());// converting decimal to whole number smaller than the number itself


			Real phi1 = grainOriNuc[indexOri + 0];
			Real PHI = grainOriNuc[indexOri + 1];
			Real phi2 = grainOriNuc[indexOri + 2];

					
				
				
						

				long DGI = -1; //Not needed for solidification

				int oi = 1; //*Index represents the affiliation to an ideal ori*/

				MPI_IO_GrainInfo grainData;// this is a struct.

				/*ix = 203;
				iy = 417;
				iz = 235;*/



				if (myNode() == cellnode) // if the generated coodinates of the nucleus lies in my process
				{
					cellP* cp = getCellPerP(ix, iy, iz);
					cellP c = *cp;
					cellP newCell = NULL;
					if (!cellExists(c))
					{
						tempCellP TCell = getTCellofCell(c);

						Real timeAtSolidification = getNucTimeInactiveCell(c);



						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]) //as: what is it doing ?
						{

							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
							// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)


							//#as in tetsing: prob.nucleation

							long len = Tcells->T_timer.size;

							for (long i = 0; i < len - 1; i++)
							{
								if ((.0000016 + timeAtSolidification) >= Tcells->T_timer.timeSeries[i] && (.0000016 + timeAtSolidification) <= Tcells->T_timer.timeSeries[i + 1])
								{
									key = i;
									break;
								}
							}
							long prevKey = 0;
							long nextKey = key + 1;

							if (key == Tcells->T_timer.size - 1) nextKey = key;

							Real ct = TCell->interpolateTemperatureSpace(ix, iy, iz, (.0000016 + timeAtSolidification), key, nextKey);  //#as in tetsing: prob.nucleation
							Real undercooling = gen.Tmelt - ct;
							//prob = .85;
							prob = this->nucprobability.calculateCumulativeProbability(undercooling, 5, .2);

							//Real prob1 = (100 * (static_cast<double>(1) / 1 * sqrt(2 * 22 / 7)) * exp(-(SQR((1375 - ct) - 5) / (2 * SQR(1))))); // gaussian

							ptest = (float)rand() / RAND_MAX;       //random number between 1 and 0

							if (prob > ptest)
							{

								if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;

								if (minNucx > ix) minNucx = ix;
								if (minNucy > iy) minNucy = iy;
								if (minNucz > iz) minNucz = iz;
								if (maxNucx < ix) maxNucx = ix;
								if (maxNucy < iy) maxNucy = iy;
								if (maxNucz < iz) maxNucz = iz;



								if (newCell)
								{
									nucleusData << 0.8 * ix << "\t" << 0.8 * iy << "\t" << 0.8 * iz << "\t" << phi1 << "\t" << PHI << "\t" << phi2 << endl; //rg.: Multiplication by 0.8 to convet to abbsolute coordiantes.
									success = 1;
									grains.Nuc++;
									prepareGrainForBroadcast(newCell, &grainData);
									

								}

							}


						}
					}
				

				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);

				if (myNode() != cellnode && success)
				{
					addGrain(&grainData);
					grains.Nuc++;
				}
			}
		} while (!success);

		bulkNucNumber--;

	}

	MPI_Barrier(MPI_COMM_WORLD);

	Real reducedTime = 0;

	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	time = reducedTime;
	Tcells->setTtimer(shortestNucTime);
	nucleusData.close(); // #rg
	write_nucleus_data_SEPI();
}

//only probablisitic nucleation 

//void caHdl::insertBulkNucleiSolidification(void)
//{
//#define NUMBER_HORIZONTAL_PASSES  5
//#define NUMBER_VERTICLE_PASSES  12//12 //for single layer it is 1
//#define SPACING_BETWEEN_HORIZONTAL_PASSES 70/0.8
//#define SPACING_BETWEEN_VERTICLE_PASSES 30/0.8
//#define LASER_START_X 50/0.8
//#define LASER_START_Y 70/0.8
//	ofstream nucleusData("NucleusData.txt", ios::out);
//	long bulkNucNumber = 23000 / 2;// 23000// number of nuclei
//	Real shortestNucTime = 1e10;//??
//	Real prob = 0;
//	Real ptest = 0;
//	long minNucx = 1000;
//	long maxNucx = 0;  //WHY max is 0 ??
//	long minNucy = 1000;
//	long maxNucy = 0;
//	long minNucz = 1000;
//	long maxNucz = 0;
//	int n = 1;
//	Real tmax = 0;
//	Real maxtime = Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]; //as:for pob nuc
//	Real tt = 0;// maxtime;//as for prob nuc
//	
//	long key;//as:added for prb nuc
//	long surfaceNucNumber = 23000 / 2;// total number
//	long bulkNucleation = 23000 / 2;
//	
//
//
//
//	
//
//
//	while (bulkNucNumber > 0) // ( TimeT < maxtime) //(bulkNucNumber > 0)
//	{
//		short success = 0;
//		do {
//			if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));
//			
//
//
//			long ix = (long)(this->r.parkMiller() * box.xPer);//parkmiller random number
//			long iz = (long)(this->r.parkMiller() * box.zPer);
//			long iy = (long)(this->r.parkMiller() * box.yPer);
//			
//
//			QUICKASSERT(iy >= 0 && iy < box.xPer);
//			QUICKASSERT(ix >= 0 && ix < box.xPer);
//			QUICKASSERT(iz >= 0 && iz < box.zPer);
//
//			int cellnode = parallelHdl->cellNode(ix, iy, iz); // find in which process that cell lies. WHICH CELL?
//
//			int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());// converting decimal to whole number smaller than the number itself
//
//
//
//
//			Real phi1 = grainOriNuc[indexOri + 0];
//			Real PHI = grainOriNuc[indexOri + 1];
//			Real phi2 = grainOriNuc[indexOri + 2];
//
//			 
//				
//
//				long DGI = -1; //Not needed for solidification
//
//				int oi = 1; //*Index represents the affiliation to an ideal ori*/
//
//				MPI_IO_GrainInfo grainData;// this is a struct.
//
//				/*ix = 203;
//				iy = 417;
//				iz = 235;*/
//
//
//
//				if (myNode() == cellnode) // if the generated coodinates of the nucleus lies in my process
//				{
//					cellP* cp = getCellPerP(ix, iy, iz);
//					cellP c = *cp;
//					cellP newCell = NULL;
//					if (!cellExists(c))
//					{
//						tempCellP TCell = getTCellofCell(c);
//
//						Real timeAtSolidification = getNucTimeInactiveCell(c);
//
//
//
//						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]) //as: what is it doing ?
//						{
//
//							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
//							// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)
//
//
//							//#as in tetsing: prob.nucleation
//
//							long len = Tcells->T_timer.size;
//
//							for (long i = 0; i < len - 1; i++)
//							{
//								if ((.0000016 + timeAtSolidification) >= Tcells->T_timer.timeSeries[i] && (.0000016 + timeAtSolidification) <= Tcells->T_timer.timeSeries[i + 1])
//								{
//									key = i;
//									break;
//								}
//							}
//							long prevKey = 0;
//							long nextKey = key + 1;
//
//							if (key == Tcells->T_timer.size - 1) nextKey = key;
//
//							Real ct = TCell->interpolateTemperatureSpace(ix, iy, iz, (.0000016 + timeAtSolidification), key, nextKey);  //#as in tetsing: prob.nucleation
//							Real undercooling = gen.Tmelt - ct;
//							//prob = .85;
//							prob = this->nucprobability.calculateCumulativeProbability(undercooling, 5, .2);  //initial undrcooling tm 5 and sigma 1
//
//							//Real prob1 = (100 * (static_cast<double>(1) / 1 * sqrt(2 * 22 / 7)) * exp(-(SQR((1375 - ct) - 5) / (2 * SQR(1))))); // gaussian
//
//							ptest = (float)rand() / RAND_MAX;       //random number between 1 and 0
//
//							if (prob > ptest)
//							{
//
//								if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;
//
//								if (minNucx > ix) minNucx = ix;
//								if (minNucy > iy) minNucy = iy;
//								if (minNucz > iz) minNucz = iz;
//								if (maxNucx < ix) maxNucx = ix;
//								if (maxNucy < iy) maxNucy = iy;
//								if (maxNucz < iz) maxNucz = iz;
//
//
//
//								if (newCell)
//								{
//									nucleusData << 0.8 * ix << "\t" << 0.8 * iy << "\t" << 0.8 * iz << "\t" << phi1 << "\t" << PHI << "\t" << phi2 << endl; //rg.: Multiplication by 0.8 to convet to abbsolute coordiantes.
//									success = 1;
//									grains.Nuc++;
//									prepareGrainForBroadcast(newCell, &grainData);
//									 
//
//								}
//
//							}
//
//
//						}
//					}
//				}
//
//				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
//				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);
//
//				if (myNode() != cellnode && success)
//				{
//					addGrain(&grainData);
//					grains.Nuc++;
//				}
//		} while (!success);
//
//			bulkNucNumber--;
//
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	Real reducedTime = 0;
//
//	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//
//	time = reducedTime;
//	Tcells->setTtimer(shortestNucTime);
//	nucleusData.close(); // #rg
//	write_nucleus_data_SEPI();
//}



/*original bulk*/

//void caHdl::insertBulkNucleiSolidification(void)
//{
//	ofstream nucleusData("NucleusData.txt", ios::out); //rg: To genearte Sepi data.
//
//	long bulkNucNumber = 23000;// 23000;
//	Real shortestNucTime = 1e10;
//
//	long minNucx = 1000;
//	long maxNucx = 0;
//	long minNucy = 1000;
//	long maxNucy = 0;
//	long minNucz = 1000;
//	long maxNucz = 0;
//
//	while (bulkNucNumber > 0)
//	{
//		short success = 0;
//		do {
//			if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));
//
//			long ix = (long)(this->r.parkMiller() * box.xPer);
//			long iy = (long)(this->r.parkMiller() * box.yPer);
//			long iz = (long)(this->r.parkMiller() * box.zPer);
//
//			QUICKASSERT(ix >= 0 && ix < box.xPer);
//			QUICKASSERT(iy >= 0 && iy < box.yPer);
//			QUICKASSERT(iz >= 0 && iz < box.zPer);
//
//			int cellnode = parallelHdl->cellNode(ix, iy, iz); // find in which process that cell lies.
//
//			int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());
//
//			Real phi1 = grainOriNuc[indexOri + 0];
//			Real PHI = grainOriNuc[indexOri + 1];
//			Real phi2 = grainOriNuc[indexOri + 2];
//
//			long DGI = -1; //Not needed for solidification
//
//			int oi = 1; //*Index represents the affiliation to an ideal ori*/
//
//			MPI_IO_GrainInfo grainData;// this is a struct.
//
//			/*ix = 250;
//			iy = 250;
//			iz = 250;*/
//
//			if (myNode() == cellnode) // if the generated coodinates of the nucleus lies in my process
//			{
//				cellP* cp = getCellPerP(ix, iy, iz);
//				cellP c = *cp;
//				cellP newCell = NULL;
//				if (!cellExists(c))
//				{
//					tempCellP TCell = getTCellofCell(c);
//
//					Real timeAtSolidification = getNucTimeInactiveCell(c);
//
//					if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1])
//					{
//						newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
//						// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)
//
//						if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;
//
//						if (minNucx > ix) minNucx = ix;
//						if (minNucy > iy) minNucy = iy;
//						if (minNucz > iz) minNucz = iz;
//						if (maxNucx < ix) maxNucx = ix;
//						if (maxNucy < iy) maxNucy = iy;
//						if (maxNucz < iz) maxNucz = iz;
//
//					}
//
//					if (newCell)
//					{
//						nucleusData << 0.8 * ix << "\t" << 0.8 * iy << "\t" << 0.8 * iz << "\t" << phi1 << "\t" << PHI << "\t" << phi2 << endl; //rg.: Multiplication by 0.8 to convet to abbsolute coordiantes.
//						success = 1;
//						grains.Nuc++;
//						prepareGrainForBroadcast(newCell, &grainData);
//					}
//				}
//			}
//
//			MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
//			MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);
//
//			if (myNode() != cellnode && success)
//			{
//				addGrain(&grainData);
//				grains.Nuc++;
//			}
//
//		} while (!success);
//
//		bulkNucNumber--;
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	Real reducedTime = 0;
//
//	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//
//	time = reducedTime;
//	Tcells->setTtimer(shortestNucTime);
//	nucleusData.close(); // #rg
//	write_nucleus_data_SEPI();
//}





//void caHdl::insertBulkNucleiSolidification(void)
//{// new code added to intoduce gradient in the microstructure.
//
//	ofstream nucleusData("NucleusData.txt", ios::out); //rg: To genearte Sepi data.
//
//	long bulkNucNumber = 23000;
//	Real shortestNucTime = 1e10;
//
//	long minNucx = 1000;
//	long maxNucx = 0;
//	long minNucy = 1000;
//	long maxNucy = 0;
//	long minNucz = 1000;
//	long maxNucz = 0;
//
//	std::vector<long>  y_values{ 70, 100 , 130, 160, 190, 220, 250, 280, 310, 340, 370, 400 };
//
//	while (bulkNucNumber > 0)
//	{
//		short success = 0;
//
//		do {
//			if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));
//
//			long ix = (long)(this->r.parkMiller() * box.xPer);
//			long iz = (long)(this->r.parkMiller() * box.zPer);
//
//			QUICKASSERT(ix >= 0 && ix < box.xPer);
//			QUICKASSERT(iz >= 0 && iz < box.zPer);
//
//			int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());
//
//			Real phi1_fixed = grainOriNuc[indexOri + 0];
//			Real PHI_fixed = grainOriNuc[indexOri + 1];
//			Real phi2_fixed = grainOriNuc[indexOri + 2];
//
//			for (int i = 0; i < y_values.size(); i++)
//			{
//				long iy = y_values[i];
//				int cellnode = parallelHdl->cellNode(ix, iy, iz); // find in which process that cell lies.
//
//				Real phi1 = phi1_fixed + (float)rand() / RAND_MAX;
//				Real PHI = PHI_fixed + (float)rand() / RAND_MAX;
//				Real phi2 = phi2_fixed + (float)rand() / RAND_MAX;
//
//				long DGI = -1; //Not needed for solidification
//
//				int oi = 1; //*Index represents the affiliation to an ideal ori*/
//
//				MPI_IO_GrainInfo grainData;// this is a struct.
//
//				if (myNode() == cellnode) // if the generated coodinates of the nucleus lies in my process
//				{
//					cellP* cp = getCellPerP(ix, iy, iz);
//					cellP c = *cp;
//					cellP newCell = NULL;
//					if (!cellExists(c))
//					{
//						tempCellP TCell = getTCellofCell(c);
//
//						Real timeAtSolidification = getNucTimeInactiveCell(c);
//
//						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1])
//						{
//							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
//						// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)
//
//							if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;
//
//							if (minNucx > ix) minNucx = ix;
//							if (minNucy > iy) minNucy = iy;
//							if (minNucz > iz) minNucz = iz;
//							if (maxNucx < ix) maxNucx = ix;
//							if (maxNucy < iy) maxNucy = iy;
//							if (maxNucz < iz) maxNucz = iz;
//
//						}
//
//						if (newCell)
//						{
//							nucleusData << 0.8*ix << "\t" << 0.8*iy << "\t" << 0.8*iz << "\t" << phi1 << "\t" << PHI << "\t" << phi2 << endl; //rg.: Multiplication by 0.8 to convet to abbsolute coordiantes.
//							success = 1;
//							bulkNucNumber--;
//							grains.Nuc++;
//							prepareGrainForBroadcast(newCell, &grainData);
//						}
//					}
//				}
//
//				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
//				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);
//
//				if (myNode() != cellnode && success)
//				{
//					addGrain(&grainData);
//					grains.Nuc++;
//				}
//			}
//
//		} while (!success);		
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	Real reducedTime = 0;
//
//	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//
//	time = reducedTime;
//	Tcells->setTtimer(shortestNucTime);
//	nucleusData.close(); // #rg
//	write_nucleus_data_SEPI();
//}


void caHdl::insertEpitaxialNucleationSolidification(void)  // line98 of cahdl.h

{
	ofstream nucleusData("EptaxNucleusData.txt", ios::out); //rg: To genearte Sepi data. making nucleas data file
#define CIRCLE 1
#define ELLIPSE 2
#define PARABOLA 3
#define CLUSTER 4
#define NUMBER_HORIZONTAL_PASSES  5
#define NUMBER_VERTICLE_PASSES  12//12 //for single layer it is 1
#define SPACING_BETWEEN_HORIZONTAL_PASSES 70/0.8
#define SPACING_BETWEEN_VERTICLE_PASSES 30/0.8
#define LASER_START_X 40/0.8
#define LASER_START_Y 70/0.8 //for multi layer
//#define LASER_START_Y 500 //considering 70 um depth for single layer
	//srand(time(0)); // to geenrate random number everytime it is run
	long bulkNucNumber = 23000;// 23000// number of nuclei

	Real prob = 0;
	Real ptest = 0;
	int n = 1;
	Real tmax = 0;
	Real maxtime = Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]; //as:for pob nuc
	Real tt = 0;// maxtime;//as for prob nuc

	int profileShape = 1;
	long surfaceNucNumber = 23000;// total number
	long fixedEpatixalNucleiNumber = 23000;

	Real shortestNucTime = 1e10;
	long numberNucleusEachLayer = fixedEpatixalNucleiNumber / NUMBER_VERTICLE_PASSES;
	long numberofNucleusEachLine = numberNucleusEachLayer / NUMBER_HORIZONTAL_PASSES;
	long minNucx = 1000;
	long maxNucx = 0;
	long minNucy = 1000;
	long maxNucy = 0;
	long minNucz = 1000;
	long maxNucz = 0;
	long key;//asaded for prb nuc
	if (profileShape == CIRCLE) {



		float radius = 60/.8;

		vector <int> y;
		vector <int> x;



		for (int i = 0; i < NUMBER_VERTICLE_PASSES; i++) {
			y.push_back(LASER_START_Y + i * SPACING_BETWEEN_VERTICLE_PASSES);
			for (int j = 0; j < NUMBER_HORIZONTAL_PASSES; j++) {
				x.push_back(LASER_START_X + j * SPACING_BETWEEN_HORIZONTAL_PASSES);

			}

		}

		while (surfaceNucNumber > 0)
		{

			short success = 0;
			do {
				if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));


				int indexY = int((fixedEpatixalNucleiNumber - surfaceNucNumber) / numberNucleusEachLayer); // this decides the maximum numbe of nucleus per layer
				if (indexY < 0) indexY = 0;
				if (indexY >= NUMBER_VERTICLE_PASSES) indexY = NUMBER_VERTICLE_PASSES - 1;
				int indexX = int((fixedEpatixalNucleiNumber - surfaceNucNumber) / numberofNucleusEachLine); // check it works ?this decides the maximum numbe of nucleus per layer
				if (indexX < 0) indexX = 0;
				if (indexX >= NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES) indexX = NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES - 1;
				QUICKASSERT(indexX <= NUMBER_VERTICLE_PASSES * NUMBER_HORIZONTAL_PASSES);


				//selction for a x y and z for inner surface
				long xo = x[indexX];
				long yo = y[indexY];

				double theta = this->r.parkMiller() * 2 * 3.14; // (rand() / (double)RAND_MAX) * 2 * 3.14;

				long ix_relative = (long)(radius * cos(theta));   // any x cordinate along the surface
				long iy_relative = abs((long)(radius * sin(theta))); // (y -y_centre)// abs used because we wan to cover only lower half

				long ix = xo + ix_relative;
				if (ix >= box.xPer)
					ix = box.xPer - 1;
				if (ix < 0)
					ix = 0; // if in case the choosen point lies outside the domain, this is due to the fact that the meltpool 
				long iy = yo - iy_relative; // because we want to cover only lower half
				if (iy >= box.yPer)
					iy = box.yPer - 1;
				long iz = (long)(this->r.parkMiller() * this->box.zPer); // (rand() / (double)RAND_MAX) * this->box.zPer;// Compute z using for entoire region


				QUICKASSERT(ix >= 0 && ix < box.xPer);
				QUICKASSERT(iy >= 0 && iy < box.yPer);
				QUICKASSERT(iz >= 0 && iz < box.zPer);

				int cellnode = parallelHdl->cellNode(ix, iy, iz);

				int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());

				Real phi1 = grainOriNuc[indexOri + 0];
				Real PHI = grainOriNuc[indexOri + 1];
				Real phi2 = grainOriNuc[indexOri + 2];

				//Real rhoTotal = 0.0;//grainOriRho[indexOri + 3];//rhoTotal could be used 

				//Real rhoGND = rhoTotal; //Not needed for solidification
				//Real rm0 = rhoTotal; //Not needed for solidification
				//Real ri0 = rhoTotal; //Not needed for solidification
				//Real rw0 = rhoTotal; //Not needed for solidification
				//int n1 = 0; //Not needed for solidification
				//int n2 = 0; //Not needed for solidification
				//int n3 = 0; //Not needed for solidification
				//Real NGLS = 0.0; //Not needed for solidification
				long DGI = -1; //Not needed for solidification

				int oi = 1; //*Index represents the affiliation to an ideal ori*/

				MPI_IO_GrainInfo grainData;

				if (myNode() == cellnode)
				{
					cellP* cp = getCellPerP(ix, iy, iz);
					cellP c = *cp;
					cellP newCell = NULL;

					if (!cellExists(c))
					{
						tempCellP TCell = getTCellofCell(c);

						Real timeAtSolidification = getNucTimeInactiveCell(c);

						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]) //as: what is it doing ?
						{

							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
							// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)


							//#as in tetsing: prob.nucleation
		
							long len = Tcells->T_timer.size;

							for (long i = 0; i < len - 1; i++)
							{
								if ((.0000016 + timeAtSolidification) >= Tcells->T_timer.timeSeries[i] && (.0000016 + timeAtSolidification) <= Tcells->T_timer.timeSeries[i + 1])
								{
									key = i;
									break;
								}
							}
							long prevKey = 0;
							long nextKey = key + 1;

							if (key == Tcells->T_timer.size - 1) nextKey = key;
							//cout << "next key is"<<nextKey;
							//long key = Tcells->T_timer.key;
							//TCellP Tcell = TCellOfActiveCell(c);
							Real ct = TCell->interpolateTemperatureSpace(ix, iy, iz, (.0000016 + timeAtSolidification), key, nextKey);  //#as in tetsing: prob.nucleation
							//if (ct > tmax) { tmax = ct; }
							Real undercooling = gen.Tmelt - ct;
							prob = this->nucprobability.calculateCumulativeProbability(undercooling,1, .2);

						//	Real prob1 = (100 * (static_cast<double>(1) / 1 * sqrt(2 * 22 / 7)) * exp(-(SQR((1375 - ct) - ) / (2 * SQR(1))))); // gaussian

							ptest = (float)rand() / RAND_MAX;       //random number between 1 and 0

							if (prob > ptest)
							{

								if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;

								if (minNucx > ix) minNucx = ix;
								if (minNucy > iy) minNucy = iy;
								if (minNucz > iz) minNucz = iz;
								if (maxNucx < ix) maxNucx = ix;
								if (maxNucy < iy) maxNucy = iy;
								if (maxNucz < iz) maxNucz = iz;



								if (newCell)
								{
									nucleusData << 0.8 * ix << "\t" << 0.8 * iy << "\t" << 0.8 * iz << "\t" << phi1 << "\t" << PHI << "\t" << phi2 << endl; //rg.: Multiplication by 0.8 to convet to abbsolute coordiantes.
									success = 1;
									grains.Nuc++;
									prepareGrainForBroadcast(newCell, &grainData);
								}

							}
							//here put here if brackets

						}
					}
				}

				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);

				if (myNode() != cellnode && success)
				{
					addGrain(&grainData);
					grains.Nuc++;
				}

			} while (!success);

			surfaceNucNumber--;
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);

	Real reducedTime = 0;

	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	time = reducedTime;
	Tcells->setTtimer(shortestNucTime);
	nucleusData.close(); // #rg
	write_nucleus_data_SEPI();


}

// Added for printing sepi data

void caHdl::write_nucleus_data_SEPI(void) {
	// write nucleus x, y, z, orientation id, theta1, theata, theata2

	ofstream nucluesDataFile("Nucleus_Data_SEPI", ios::out);

	nucluesDataFile << "Nucleus No." << "\t" << "x_coor" << "\t" << "y_corr" << "\t" << "z_corr" << "\t" << "Orientation_id" << "\t" << "phi_1" << "\t" << "PHI" << "\t" << "phi_2" << "\t" << endl;

	if (!nucluesDataFile) {
		cout << "File could not be created" << endl;
	}
	for (int i = 0; i < this->grains.Nuc + 1; i++) {//+1 for undeformed grain powder

		long dummy_orientation_index = grains.all[i]->oriIndex;

		float local_xcorr = this->grains.all[i]->mySpace.xleft;
		float local_ycorr = this->grains.all[i]->mySpace.yleft;
		float local_zcorr = this->grains.all[i]->mySpace.zleft;

		double local_phi_1 = this->grains.orientations[dummy_orientation_index]->phi1;
		double local_phi = this->grains.orientations[dummy_orientation_index]->PHI;
		double local_phi_2 = this->grains.orientations[dummy_orientation_index]->phi2;

		//nucluesDataFile << i << "\t" << local_xcorr << "\t" << local_ycorr << "\t" << local_zcorr << "\t"<< dummy_orientation_index << "\t" << local_phi_1 << "\t" << local_phi << "\t" << local_phi_2 << endl;
		nucluesDataFile << i << "\t" << 0.8 * local_xcorr << "\t" << 0.8 * local_ycorr << "\t" << 0.8 * local_zcorr << "\t" << dummy_orientation_index << "\t" << local_phi_1 << "\t" << local_phi << "\t" << local_phi_2 << endl;

	}

	nucluesDataFile.close();
}

// Till here upper function was added to print sepi data.

cellP* caHdl::getCellPerP(long ix, long iy, long iz) // gives the local id of the cell with in the process.
{
	QUICKASSERT(myNode() == parallelHdl->cellNode(ix, iy, iz));
	return &cells.all[(iz - z0) * cells.data.zFac + (iy - y0) * cells.data.xPer + (ix - x0)]; //# as: what is this ?
}

cellP caHdl::newNucleusSol(Real p1, Real ps, Real p2, long ix, long iy, long iz, TCellP Tcell, long DeformGrainIndex, int Oriindex, Real TimeToMelt)
{
	cell pseudoNb;
	pseudoNb.ix = (short)ix;
	pseudoNb.iy = (short)iy;
	pseudoNb.iz = (short)iz;
	pseudoNb.dataHolder = 0;
	pseudoNb.rxFrac = 0;

	cellP nucCell = infectCell(&pseudoNb, 0, 0, 0, 26, 0.0);	//Here the cell is created but the rxGrain and P assume false values

	if (nucCell == NOTNEWCELLCREATED) return NULL;

	grainP nucGrain = new grain(&(this->grains), p1, ps, p2, 0, 0, 0, 0, 0, 0, 0, Oriindex, 1, DeformGrainIndex);//Check: 5th to 11th argument not needed                                                                                                                                                                                                                                                                                                                                                          (&(this->grains), p1, ps, p2, 0, 0, 0, 0, 0, 0, 0, Oriindex, 1, DeformGrainIndex);//Check: 5th to 11th argument not needed
	nucGrain->rhoTotal = 0; //Check: It could be used

	nucGrain->mySpace.xleft = nucGrain->mySpace.xright = (short)ix;
	nucGrain->mySpace.yleft = nucGrain->mySpace.yright = (short)iy;
	nucGrain->mySpace.zleft = nucGrain->mySpace.zright = (short)iz;

	nucCell->dataHolder &= 0xFFFFF8000000001F; //###LB: Mask to avoid false values

	uint_fast64_t gpos = (uint_fast64_t)nucGrain->arrayPos;
	gpos <<= 22; // storing rx grain
	nucCell->dataHolder |= gpos;

	uint_fast64_t TCpos = (uint_fast64_t)Tcell->id;
	TCpos <<= 43; // storing temp cell.
	nucCell->dataHolder |= TCpos;

	grainP rxGrain = grainOfCell(nucCell);

	//Check: Instead of P encode time to melt

	uint_fast64_t timeDiscrete = 0xFFFF;

	if ((long)TimeToMelt != LIQUID && (long)TimeToMelt != NO_MELT) //as:added
		timeDiscrete = (uint_fast64_t)(TimeToMelt * Tcells->timeEncodeFactor);

	timeDiscrete <<= 5;// storing nucleation time.

	nucCell->dataHolder |= timeDiscrete;

	QUICKASSERT(rxGrain == nucGrain);

	QUICKASSERT(nucCell);

	nucGrain->CellCount = 0;
	nucGrain->initialCellNumber = 0;

	Real testTime = Tcells->getCellNucleationTime(nucCell->dataHolder);

	return nucCell;
}

void caHdl::prepareGrainForBroadcast(cellP newCell, MPI_IO_GrainInfo* grainData)
{
	grainP g = grainOfCell(newCell);		//The grain that is broadcasted is the new grain, deformed grains are already in the lists
	if (!g) err->reportError(ERRTXT("Grain must exist"));

	grainData->alreadyNuc = g->alreadyNuc;
	grainData->CellCount = g->CellCount;
	grainData->initialCellNumber = g->initialCellNumber;
	if (grainData->initialCellNumber < 0) err->reportError(ERRTXT("Must be positive or zero"));
	grainData->colB = g->colB;
	grainData->colG = g->colG;
	grainData->colR = g->colR;
	grainData->DeformGrainIndex = g->DeformGrainIndex;
	grainData->NGLS = g->NGLS;
	grainData->nuc1 = g->nuc1;
	grainData->nuc2 = g->nuc2;
	grainData->nuc3 = g->nuc3;
	grainData->Oriindex = g->Oriindex;
	grainData->oriIndex = g->oriIndex;

	/*grainData->phi1_TB = g->phi1_TB;
	grainData->phi2_TB = g->phi2_TB;
	grainData->PHI_TB = g->PHI_TB;*/

	grainData->rhoGND = g->rhoGND;
	grainData->rhoI = g->rhoI;
	grainData->rhoIntern0 = g->rhoIntern0;
	grainData->rhoM = g->rhoM;
	grainData->rhoMobil0 = g->rhoMobil0;
	grainData->rhoTotal = g->rhoTotal;
	grainData->rhoW = g->rhoW;
	grainData->rhoWall0 = g->rhoWall0;
	grainData->RVFac3 = g->RVFac3;
	grainData->scatter = g->scatter;
	grainData->phi1 = grains.orientations[g->oriIndex]->phi1;
	grainData->PHI = grains.orientations[g->oriIndex]->PHI;
	grainData->phi2 = grains.orientations[g->oriIndex]->phi2;

	grainData->arrayPos = g->arrayPos;
	grainData->xleft = g->mySpace.xleft;
	grainData->xright = g->mySpace.xright;
	grainData->yleft = g->mySpace.yleft;
	grainData->yright = g->mySpace.yright;
	grainData->zleft = g->mySpace.zleft;
	grainData->zright = g->mySpace.zright;
}

void caHdl::prepareGrainForBroadcast(long grainsArrayPos, MPI_IO_GrainInfo* grainData)
{
	grainP g = grains.all[grainsArrayPos];		//The grain that is broadcasted is the new grain, deformed grains are already in the lists
	if (!g)
		err->reportError(ERRTXT("Grain must exist"));

	grainData->alreadyNuc = g->alreadyNuc;
	grainData->CellCount = g->CellCount;
	grainData->initialCellNumber = g->initialCellNumber;

	if (grainData->initialCellNumber < 0)
		err->reportError(ERRTXT("Must be positive or zero"));

	grainData->colB = g->colB;
	grainData->colG = g->colG;
	grainData->colR = g->colR;
	grainData->DeformGrainIndex = g->DeformGrainIndex;
	grainData->NGLS = g->NGLS;
	grainData->nuc1 = g->nuc1;
	grainData->nuc2 = g->nuc2;
	grainData->nuc3 = g->nuc3;
	grainData->Oriindex = g->Oriindex;
	grainData->oriIndex = g->oriIndex;
	/*grainData->phi1_TB = g->phi1_TB;
	grainData->phi2_TB = g->phi2_TB;
	grainData->PHI_TB = g->PHI_TB;*/
	grainData->rhoGND = g->rhoGND;
	grainData->rhoI = g->rhoI;
	grainData->rhoIntern0 = g->rhoIntern0;
	grainData->rhoM = g->rhoM;
	grainData->rhoMobil0 = g->rhoMobil0;
	grainData->rhoTotal = g->rhoTotal;
	grainData->rhoW = g->rhoW;
	grainData->rhoWall0 = g->rhoWall0;
	grainData->RVFac3 = g->RVFac3;
	grainData->scatter = g->scatter;
	grainData->phi1 = grains.orientations[g->oriIndex]->phi1;
	grainData->PHI = grains.orientations[g->oriIndex]->PHI;
	grainData->phi2 = grains.orientations[g->oriIndex]->phi2;

	grainData->arrayPos = g->arrayPos;
	grainData->xleft = g->mySpace.xleft;
	grainData->xright = g->mySpace.xright;
	grainData->yleft = g->mySpace.yleft;
	grainData->yright = g->mySpace.yright;
	grainData->zleft = g->mySpace.zleft;
	grainData->zright = g->mySpace.zright;
}

//###from performance consideration additional if but is sufficient for the time being as initialization phase is short and will sooner or later replaced by a new microstructure generator/input framework
//###should be avoided during setting millions of nuclei

grainP caHdl::addGrain(MPI_IO_GrainInfo* g)
{
	grainP nucGrain = new grain(&(this->grains), g->phi1, g->PHI, g->phi2, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, -1);

	nucGrain->alreadyNuc = g->alreadyNuc;
	nucGrain->CellCount = g->CellCount;
	nucGrain->initialCellNumber = g->initialCellNumber;
	if (nucGrain->initialCellNumber < 0) err->reportError(ERRTXT("Must be positive or zero"));
	nucGrain->colB = g->colB;
	nucGrain->colG = g->colG;
	nucGrain->colR = g->colR;
	nucGrain->DeformGrainIndex = g->DeformGrainIndex;
	nucGrain->NGLS = g->NGLS;
	nucGrain->nuc1 = g->nuc1;
	nucGrain->nuc2 = g->nuc2;
	nucGrain->nuc3 = g->nuc3;
	nucGrain->Oriindex = g->Oriindex;

	if (nucGrain->oriIndex != g->oriIndex)
		err->reportError(ERRTXT("Same grain in different processors must point to same ori"));
	if (nucGrain->arrayPos != g->arrayPos)
		err->reportError(ERRTXT("Position in array must be the same"));

	/*nucGrain->phi1_TB = g->phi1_TB;
	nucGrain->phi2_TB = g->phi2_TB;
	nucGrain->PHI_TB = g->PHI_TB;*/

	nucGrain->rhoGND = g->rhoGND;
	nucGrain->rhoI = g->rhoI;
	nucGrain->rhoIntern0 = g->rhoIntern0;
	nucGrain->rhoM = g->rhoM;
	nucGrain->rhoMobil0 = g->rhoMobil0;
	nucGrain->rhoTotal = g->rhoTotal;
	nucGrain->rhoW = g->rhoW;
	nucGrain->rhoWall0 = g->rhoWall0;
	nucGrain->RVFac3 = g->RVFac3;
	nucGrain->scatter = g->scatter;

	nucGrain->mySpace.xleft = g->xleft;
	nucGrain->mySpace.xright = g->xright;
	nucGrain->mySpace.yleft = g->yleft;
	nucGrain->mySpace.yright = g->yright;
	nucGrain->mySpace.zleft = g->zleft;
	nucGrain->mySpace.zright = g->zright;
	return nucGrain;
}

void caHdl::insertSurfaceNucleiSolidification(void)
{
#define CIRCLE 1
#define ELLIPSE 2
#define PARABOLA 3
#define CLUSTER 4

#define NUMBER_HORIZONTAL_PASSES  5
#define NUMBER_VERTICLE_PASSES  12

#define SPACING_BETWEEN_HORIZONTAL_PASSES 100/0.8
#define SPACING_BETWEEN_VERTICLE_PASSES 30/0.8

#define LASER_START_X 50/0.8
#define LASER_START_Y 70/0.8

	int profileShape = 4;
	long surfaceNucNumber = 23000;// total number
	long fixedSurfaceNucleiNumber = 23000;

	Real shortestNucTime = 1e10;

	long minNucx = 1000;
	long maxNucx = 0;
	long minNucy = 1000;
	long maxNucy = 0;
	long minNucz = 1000;
	long maxNucz = 0;

	if (profileShape == CIRCLE) {

		float x;
		float y;
		float z;

		float radius = 50;

		int x_centre[NUMBER_HORIZONTAL_PASSES] = { 0 };
		int y_centre[NUMBER_VERTICLE_PASSES] = { 0 };

		for (int i = 0; i < NUMBER_HORIZONTAL_PASSES; i++) {
			x_centre[i] = LASER_START_X + i * SPACING_BETWEEN_HORIZONTAL_PASSES;
		}

		for (int j = 0; j < NUMBER_VERTICLE_PASSES; j++) {
			y_centre[j] = LASER_START_Y + j * SPACING_BETWEEN_VERTICLE_PASSES;
		}

		while (surfaceNucNumber > 0)
		{
			short success = 0;
			do {
				if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));

				int dummy_x_centre_index = (int)(this->r.parkMiller() * NUMBER_HORIZONTAL_PASSES); // choose any x centre
				int dummy_y_centre_index = (int)(this->r.parkMiller() * NUMBER_VERTICLE_PASSES); // choose any y centre

				long xo = x_centre[dummy_x_centre_index];
				long yo = y_centre[dummy_y_centre_index];

				double theta = this->r.parkMiller() * 2 * 3.14;

				long ix_relative = (long)(radius * cos(theta));   // (x - x_centre))
				long iy_relative = abs((long)(radius * sin(theta))); // (y -y_centre)// abs used because we wan to cover oknly lower half

				long ix = xo + ix_relative;
				long iy = yo - iy_relative; // because we want to cover only lower half
				long iz = (long)(this->r.parkMiller() * this->box.zPer);

				QUICKASSERT(ix >= 0 && ix < box.xPer);
				QUICKASSERT(iy >= 0 && iy < box.yPer);
				QUICKASSERT(iz >= 0 && iz < box.zPer);

				int cellnode = parallelHdl->cellNode(ix, iy, iz);

				int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());

				Real phi1 = grainOriNuc[indexOri + 0];
				Real PHI = grainOriNuc[indexOri + 1];
				Real phi2 = grainOriNuc[indexOri + 2];

				//Real rhoTotal = 0.0;//grainOriRho[indexOri + 3];//rhoTotal could be used 

				//Real rhoGND = rhoTotal; //Not needed for solidification
				//Real rm0 = rhoTotal; //Not needed for solidification
				//Real ri0 = rhoTotal; //Not needed for solidification
				//Real rw0 = rhoTotal; //Not needed for solidification
				//int n1 = 0; //Not needed for solidification
				//int n2 = 0; //Not needed for solidification
				//int n3 = 0; //Not needed for solidification
				//Real NGLS = 0.0; //Not needed for solidification
				long DGI = -1; //Not needed for solidification

				int oi = 1; //*Index represents the affiliation to an ideal ori*/

				MPI_IO_GrainInfo grainData;

				if (myNode() == cellnode)
				{
					cellP* cp = getCellPerP(ix, iy, iz);
					cellP c = *cp;
					cellP newCell = NULL;

					if (!cellExists(c))
					{
						tempCellP TCell = getTCellofCell(c);

						Real timeAtSolidification = getNucTimeInactiveCell(c);

						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1])
						{
							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
							// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)

							if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;

							if (minNucx > ix) minNucx = ix;
							if (minNucy > iy) minNucy = iy;
							if (minNucz > iz) minNucz = iz;
							if (maxNucx < ix) maxNucx = ix;
							if (maxNucy < iy) maxNucy = iy;
							if (maxNucz < iz) maxNucz = iz;

						}

						if (newCell)
						{
							success = 1;
							grains.Nuc++;
							prepareGrainForBroadcast(newCell, &grainData);
						}
					}
				}

				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);

				if (myNode() != cellnode && success)
				{
					addGrain(&grainData);
					grains.Nuc++;
				}

			} while (!success);

			surfaceNucNumber--;
		}
	}

	if (profileShape == ELLIPSE) {}

	if (profileShape == PARABOLA) {}

	if (profileShape == CLUSTER)
	{
		std::vector<int> y_values;
		for (auto i = 0; i < NUMBER_VERTICLE_PASSES; i++)
		{
			y_values.push_back(LASER_START_Y + i * SPACING_BETWEEN_VERTICLE_PASSES);
			if (y_values[i] < 0) y_values[i] = 0;
			if (y_values[i] >= box.yPer) y_values[i] = box.yPer - 1;
		}

		long numberNucleusEachLayer = fixedSurfaceNucleiNumber / NUMBER_VERTICLE_PASSES;

		while (surfaceNucNumber > 0)
		{
			short success = 0;
			do {
				if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));

				int index = int((fixedSurfaceNucleiNumber - surfaceNucNumber) / numberNucleusEachLayer);
				if (index < 0) index = 0;
				if (index >= NUMBER_VERTICLE_PASSES) index = NUMBER_VERTICLE_PASSES - 1;

				long iy = y_values[index];
				long ix = (long)(this->r.parkMiller() * box.xPer);
				long iz = (long)(this->r.parkMiller() * box.zPer);

				QUICKASSERT(ix >= 0 && ix < box.xPer);
				QUICKASSERT(iy >= 0 && iy < box.yPer);
				QUICKASSERT(iz >= 0 && iz < box.zPer);

				int cellnode = parallelHdl->cellNode(ix, iy, iz);

				int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());

				Real phi1 = grainOriNuc[indexOri + 0];
				Real PHI = grainOriNuc[indexOri + 1];
				Real phi2 = grainOriNuc[indexOri + 2];

				//Real rhoTotal = 0.0;//grainOriRho[indexOri + 3];//rhoTotal could be used 

				//Real rhoGND = rhoTotal; //Not needed for solidification
				//Real rm0 = rhoTotal; //Not needed for solidification
				//Real ri0 = rhoTotal; //Not needed for solidification
				//Real rw0 = rhoTotal; //Not needed for solidification
				//int n1 = 0; //Not needed for solidification
				//int n2 = 0; //Not needed for solidification
				//int n3 = 0; //Not needed for solidification
				//Real NGLS = 0.0; //Not needed for solidification
				long DGI = -1; //Not needed for solidification

				int oi = 1; //*Index represents the affiliation to an ideal ori*/

				MPI_IO_GrainInfo grainData;

				if (myNode() == cellnode)
				{
					cellP* cp = getCellPerP(ix, iy, iz);
					cellP c = *cp;
					cellP newCell = NULL;

					if (!cellExists(c))
					{
						tempCellP TCell = getTCellofCell(c);

						Real timeAtSolidification = getNucTimeInactiveCell(c);

						if (timeAtSolidification < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1])
						{
							newCell = newNucleusSol(phi1, PHI, phi2, ix, iy, iz, TCell, DGI, oi, timeAtSolidification); //All this arguments are not needed.
							// phi1 = Real p1 || PHI = Real ps || phi2 = Real p2 || ix= long ix || iy = long iy ||iz = long iz|| TCell = TCellP Tcell ||DGI = long DeformGrainIndex || oi = int Oriindex ||  timeAtSolidification = Real TimeToMelt)

							if (timeAtSolidification < shortestNucTime) shortestNucTime = timeAtSolidification;

							if (minNucx > ix) minNucx = ix;
							if (minNucy > iy) minNucy = iy;
							if (minNucz > iz) minNucz = iz;
							if (maxNucx < ix) maxNucx = ix;
							if (maxNucy < iy) maxNucy = iy;
							if (maxNucz < iz) maxNucz = iz;

						}

						if (newCell)
						{
							success = 1;
							grains.Nuc++;
							prepareGrainForBroadcast(newCell, &grainData);
						}
					}
				}

				MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
				MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);

				if (myNode() != cellnode && success)
				{
					addGrain(&grainData);
					grains.Nuc++;
				}

			} while (!success);

			surfaceNucNumber--;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	Real reducedTime = 0;

	MPI_Allreduce(&shortestNucTime, &reducedTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	time = reducedTime;
	Tcells->setTtimer(shortestNucTime);
	//nucleusData.close(); // #rg
	//write_nucleus_data_SEPI();

}

//void caHdl::insertSurfaceNucleiSolidification(void)
//{
//	long surfaceNucNumber = 10000;
//	//We will consider to begin with the positive z direction as being an open surface
//
//	while (surfaceNucNumber > 0)
//	{
//		short success = 0;
//		do {
//			if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));
//
//			long ix = (long)(this->r.parkMiller() * box.xPer);
//			long iy = (long)(this->r.parkMiller() * box.yPer);
//			long iz = (long)(this->r.parkMiller() * box.zPer);
//
//			Real coin = r.parkMiller();
//
//			if (coin >= 0 && coin <= 0.2) ix = 0;				//This decides on which surface to nucleate
//			else if (coin > 0.2 && coin <= 0.4) ix = box.xPer - 1;
//			else if (coin > 0.4 && coin <= 0.6) iy = 0;
//			else if (coin > 0.6 && coin <= 0.8) iy = box.yPer - 1;
//			else iz = 0;
//
//			QUICKASSERT(ix >= 0 && ix < box.xPer);
//			QUICKASSERT(iy >= 0 && iy < box.yPer);
//			QUICKASSERT(iz >= 0 && iz < box.zPer);
//
//			int cellnode = parallelHdl->cellNode(ix, iy, iz);
//
//			int indexOri = 3 * floor(numberInputNucOris * myRandom.parkMiller());
//
//			Real phi1 = grainOriNuc[indexOri + 0];
//			Real PHI = grainOriNuc[indexOri + 1];
//			Real phi2 = grainOriNuc[indexOri + 2];
//
//			Real rhoTotal = 0.0;//grainOriRho[indexOri + 3];//rhoTotal could be used 
//
//			Real rhoGND = rhoTotal; //Not needed for solidification
//			Real rm0 = rhoTotal; //Not needed for solidification
//			Real ri0 = rhoTotal; //Not needed for solidification
//			Real rw0 = rhoTotal; //Not needed for solidification
//			int n1 = 0; //Not needed for solidification
//			int n2 = 0; //Not needed for solidification
//			int n3 = 0; //Not needed for solidification
//			Real NGLS = 0.0; //Not needed for solidification
//			long DGI = -1; //Not needed for solidification
//
//			int oi = 1; //orientation index
//
//			MPI_IO_GrainInfo grainData;
//
//			if (myNode() == cellnode)
//			{
//				cellP c = getCellPer(ix, iy, iz);
//				cellP newCell = NULL;
//
//				if (!c)
//				{
//					Real Tc = (Real)cellTemperature(*getCellPerP(ix, iy, iz)) * TempFactor;
//					newCell = newNucleusMG(phi1, PHI, phi2, ix, iy, iz, Tc, rhoTotal, rhoGND, rm0, ri0, rw0, n1, n2, n3, NGLS, DGI, oi);
//					if (newCell)
//					{
//						success = 1;
//						grains.Nuc++;
//						prepareGrainForBroadcast(newCell, &grainData);
//					}
//				}
//			}
//
//			MPI_Bcast(&success, 1, MPI_INT, cellnode, MPI_COMM_WORLD);
//			MPI_Bcast(&grainData, 1, MPI_IO_GrainInfo_Type, cellnode, MPI_COMM_WORLD);
//
//			if (myNode() != cellnode && success)
//			{
//				addGrain(&grainData);
//				grains.Nuc++;
//			}
//
//		} while (!success);
//
//		surfaceNucNumber--;
//	}
//}

cellP caHdl::getCellPer(long ix, long iy, long iz)
{
	QUICKASSERT(myNode() == parallelHdl->cellNode(ix, iy, iz));// get in which process that ca cell is present
	cellP c = *getCellPerP(ix, iy, iz);
	if (cellExists(c))		return c;
	return NULL;
}

cellP caHdl::newNucleusMG(Real p1, Real ps, Real p2, long ix, long iy, long iz, Real Temperature, Real rT, Real rGND, double rm0, double ri0, double rw0, int n1, int n2, int n3, Real NGLS, long DeformGrainIndex, int Oriindex)
{ // #as:this function is being called in surfacenucleation part 
	cell pseudoNb;
	pseudoNb.ix = (short)ix;
	pseudoNb.iy = (short)iy;
	pseudoNb.iz = (short)iz;
	pseudoNb.dataHolder = 0;

	cellP nucCell = infectCell(&pseudoNb, 0, 0, 0, 26, 0.0);	//Here the cell is created but the rxGrain and P assume false values
	//(cellP infector, long dx, long dy, long dz, int infectorvictimindex, Real overshoot)
	if (nucCell == NOTNEWCELLCREATED) return NULL;

	grainP nucGrain = new grain(&(this->grains), p1, ps, p2, rGND, rm0, ri0, rw0, n1, n2, n3, Oriindex, 1, DeformGrainIndex);
	nucGrain->rhoTotal = rT;

	nucGrain->mySpace.xleft = nucGrain->mySpace.xright = (short)ix;
	nucGrain->mySpace.yleft = nucGrain->mySpace.yright = (short)iy;
	nucGrain->mySpace.zleft = nucGrain->mySpace.zright = (short)iz;

	nucCell->dataHolder &= 0xFFFFF8000000001F; //###LB: Mask to avoid false values

	nucCell->T = Temperature;

	uint_fast64_t gpos = (uint_fast64_t)nucGrain->arrayPos;
	gpos <<= 22;
	nucCell->dataHolder |= gpos;

	grainP rxGrain = grainOfCell(nucCell);

	double P = 1.0;

	uint_fast64_t Pbit = 0;

	if (P < 0.0)
		Pbit = 0x200000;	//sign bit is assigned
	else
		Pbit = ((uint_fast64_t)(P * (_MOBFACUNIT))) << 5;

	nucCell->dataHolder |= Pbit;

	QUICKASSERT(rxGrain == nucGrain);

	QUICKASSERT(nucCell);

	nucGrain->CellCount = 0;
	nucGrain->initialCellNumber = 0;

	return nucCell;
}

void caHdl::determineInitialTexture(void)
{
	for (int iz = z0; iz < z0 + zPer; iz++) //  for all the ca cells in my process
		for (int iy = y0; iy < y0 + yPer; iy++)
			for (int ix = x0; ix < x0 + xPer; ix++)
			{

				cellP c = getCellPer(ix, iy, iz);  //Check this because it should not work. It works only by coincidence.
				grainP g = NULL;

				if (c) g = grainOfCell(c);
				else	g = getOrigGrain(c);

				int oriIndex = g->oriIndex;
				grains.orientations[oriIndex]->count++;
			}
}

void caHdl::startSimulation(void)
{
	long noActionCount = 0;
	long stillRXing = 1;

	char 
		active = _YES_;
	char sendFlag = NOTHING;

	deltaTime = initialDeltaTime(); //Check: This function must be modified to account for the physics of solidification
	// write_voxeldata_coloring_ipfz();
	// draw_pole_figure();
	do
	{
		if (noActionCount < NOACTIONOUTCOUNT) //TIMESLOTCOUNT ) NOACTIONOUTCOUNT is 100000
		{
			if (grains.Nuc > 0.0)
			{
				stillRXing = calcGrowthStepMG();

				if (stillRXing)		noActionCount = 0;
				else					noActionCount++;

				updateFullRX();
			}
		}

		if (noActionCount >= NOACTIONOUTCOUNT) //TIMESLOTCOUNT ) 
			active = _NO_;

		char activeOp = _NO_;
		MPI_Allreduce(&active, &activeOp, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);
		active = activeOp;

		MPI_Allreduce(&(this->somethingToSend), &sendFlag, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);

		if (active)
			noActionCount = 0;

		if (transfer.now() && sendFlag)
			sendReceive();

		if (somethingReceived) {
			distributeReceivedCells();
			noActionCount = 0;
		}

		if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) err->reportError(ERRTXT("MPI Barrier was not successful"));

		Output();
		time += gen.delta;
	//time += deltaTime;
		
		//time += gen.delta;
		
	/*	if (parm.vMax <= 0) {   //as: aadded to consider maximum crystal growth
			time += deltaTime;
		// //time += (.2 * cells.data.size / parm.vMax);
			
		//	//cout << "vmax"<<parm.vMax;
		}
		else
		{
			time += gen.delta;
			//time += deltaTime;
			

		//time += MIN (deltaTime,(cells.data.size / parm.vMax));//gen.delta; //(.2 * cells.data.size / parm.vMax);
			//cout << "new delta time" << (.2 * cells.data.size / parm.vMax) << endl; //debug
			//cout << "check time" << time << endl;
		} // debug
		*/

		if (time > Tcells->T_timer.timeSeries[Tcells->T_timer.key + 1]) //Check: Condition for last value
			Tcells->T_timer.key++;

	} while (active && time < Tcells->T_timer.timeSeries[Tcells->T_timer.size - 1]); //Check: Include conditions when tíme in CA is equal or larger than max time from FEM

	Real currentTime = this->time;
	Real maxTime = 0.0;

	MPI_Allreduce(&currentTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	cout << "writing porosity positions" << endl;
	write_prorosity();
	cout << "writing voxel datd" << endl;
	//write_voxeldata_oris(); //LB2019: For exporting the oris as voxel data.
	write_voxeldata_coloring_ipfz();
	cout << "writing voxel orientation" << endl;

	write_xsection_coloring_ipfz(); 
	orientations_SEPI();
	writingVoxelOrientationID_SEPI();
	write_final_state_SEPI();
	draw_pole_figure();

	this->time = maxTime;
	int ng = 0;
	long totalOfCells = 0;
}

Real caHdl::initialDeltaTime(void)
{
	/*if (gen.delta <=0 || gen.delta >gen.ASPIREDFILLPERSTEP * cells.data.size)
		return(gen.ASPIREDFILLPERSTEP * cells.data.size, gen.delta);
	else*/
//return MIN (gen.ASPIREDFILLPERSTEP * cells.data.size, gen.delta); //gen.delta; // //.2*.000008=.00000016 refer line 531
	return gen.ASPIREDFILLPERSTEP * cells.data.size;
}

long caHdl::calcGrowthStepMG(void) //Changed
{
	cellP c; // create a new cell struct
	cellP prev = NULL; // create another cell struct.

	Real Tm = gen.Tmelt; // #as: Melting temperature; it should be read from the input file it is in general data
	Real eta1 = gen.eta1; //1.49e-6; // m/(sK^3) -> It must be read from input file //1.49e-6 #as reeading from input file
	Real eta2 = gen.eta2; //2.9e-6; // m/(sK^2) -> It must be read from input file //2.9e-6#as reeading from input file 
	//Real dTc = gen.Uc; // #as: reading from input file
	Real dTc = 0;
	Real dTm = gen.Tcrit_mean;
	Real dTe = gen.Tcrit_sigma;
	Real totalRXsum = (Real)cells.data.retiredCount;
	Real testDt = gen.ASPIREDFILLPERSTEP * cells.data.size; //.0008*.2 cells.data.size is .8 micron
	//Real _v0 = .2 * cells.data.size / parm.vMax; //#as added for function check
	 Real _v0 = testDt / cells.data.size;  // this is aspiredfilledperstep i.e .2
	Real V = 0;
	//gen.delta = .2 * .0004 / 500;


	//added for checking this function as
 //   ofstream growth;//added for checking this function as
	//growth.open("calcGrowthStepMG.txt");//added for checking this function as
	//if (!growth) cout << "error line 2313" << endl;//added for checking this function as
	//added for checking this function as


	//set geometrical factors fgeo[infectorvictimindex]

	long somethingHasMoved = 0;

	Real maxCellLoopTime = 1e30; //as: what is this ?

	//added for checking this function as
	//growth << "testDt is" << testDt << "\t" << "_vo is" << _v0 << "totalRXsum is"<< totalRXsum<< endl;//added for checking this function as
	//growth << "c \t\t\t cellNucTime \t\t\t temperature \t\t\t speed \t\t\t jump \t\t\t vpair \t\t\ rxFrac" << endl;//added for checking this function as

	//added for checking this function as
	for (c = cells.first[0]; c; c = cells.next) //Kette aller aktuell aktiven Zellen, die erste Zelle aus der Geschwindigkeitsklasse, dann ist dynamisch die nächste Zelle in der Hauptschleife die nächste Zelle auf die die aktuellen Zelle zeigt, solange bis das ein NULL ist
	{	//Chain of all currently active cells, the first cell from the speed class, then dynamically the next cell in the main loop is the next cell to which the current cell points until this is a NULL
		cells.next = c->next;
		//cout << c->ix << "\t" << c->iy << "\t" << c->iz << endl;
		Real T = 0;
		Real dT = 0;
		Real p = 0;

		Real cellNucTime = Tcells->getCellNucleationTime(c->dataHolder);
		//Real testTemp = Tcells->getCellNucleationTime(c);//testTemp never used.
		//cout << cellNucTime << "\t" << testTemp<< endl;
		
		//growth << c << "\t\t\t" << cellNucTime << "\t\t\t";//added for checking this function as

		if (maxCellLoopTime > cellNucTime) maxCellLoopTime = cellNucTime;

		if (time >= cellNucTime)
		{
 			T = currentCellTemperature(c, this->time, Tcells->T_timer.key);
			dT = Tm - T;
			if (dT > 0)
			{
				dTc = 1*(1 / 1*sqrt(2 * 22 / 7)) * exp(-(SQR(dT - dTm) / (2 * SQR(dTe)))); // added for critical temperature
				if (dTc >1) 
				//cout << "dtc" << dTc << endl;
				Real prob = (float)rand()/RAND_MAX;//random number between 1 and 0
				if (dT > dTc)
				{
					p = eta2 * SQR(dT) + eta1 * CUBE(dT); // the velocity of 100 direction
					//V = p;
					parm.vMax = p;
					if (V < p) {
						V = p;
						  parm.vMax= V;
					}
					
					//cout << "timedelta\t" << gen.delta << endl;
					
				}

				/*if (dT < dTc)
					cout << dT << "\t" << "dtc:: " << dTc << endl;*/
			}

			//if (dT >= dTc)  //#original it is there
			//{
			//	p = eta2 * SQR(dT) + eta1 * CUBE(dT);// the velocity of 100 direction //as: for fixed crit temp 9.5
			//	//growth << p << "\t\t\t"; //added for checking this function as
			//}
			//cout << T << "\t" << endl;
			////////////#as debugging
			//cout << c << "\t\t check in line no 2361" << endl;
			//cout << c->ix << "\t\t" << c->iy << "\t\t" << c->iz<<"\t\t check in line no 2361" << endl;

			////////////#as debugging
			//growth << T << "\t\t\t"; //added for checking this function as
		}

		Real vpair = 0.0;//as what is this ?
		Real jump = TIMESLOTCOUNT;

		//growth << jump << "\t\t\t"; //added for checking this function as


		QUICKASSERT(c->rxFrac >= 0);

		QUICKASSERT(p >= 0);

		if (p > 0)
		{
			int infector = cellInfector(c);

			vpair = parm.fgeo[infector] * p;  // p is already a velocity
			 
			//growth << vpair << "\t\t\t"; //added for checking this function as9aspire

			//double dRX = vpair * _v0 + DRX; //original // _v0=testDt/cekks.data.size i.e _v0=.2 , aspiredfilledperstep i.e .2; vpair is velocity*fgeo factor; ///DOUBT _v0 should be deltaT
			
			double dRX =  (vpair * gen.delta/ cells.data.size) + DRX; //#as:added to control the infection of the cell
			
			//Real dtFac = gen.ASPIREDFILLPERSTEP / dRX; //original
			Real dtFac = (vpair * gen.delta / cells.data.size) / dRX;
			//Real dtFac =  cells.data.size/ dRX; ///(vpair * gen.delta/ cells.data.size) / dRX; //#as new check
			//Real dtFac = (vpair * gen.delta/ cells.data.size) / dRX;//#as new check
			if (dRX > 0.0000001)      somethingHasMoved = 1;

			jump = (long)(dtFac + 0.5);		// long rundet ab
			if (jump < 1) jump = 1;

			if (jump > TIMESLOTCOUNT)	jump = TIMESLOTCOUNT;


			c->rxFrac += dRX * jump;  //resembles dRX = v*deltat/d

			//growth << c->rxFrac ; //added for checking this function as
		}

		totalRXsum += c->rxFrac;

		if (c->rxFrac < 1.0) //if the cell is not fully recrystallized
		{
			long j = (long)jump;
			c->next = cells.first[j];
			cells.first[j] = c;
			if (prev) {
				prev->next = cells.next;
			}
			else {
				cells.first[0] = cells.next;
			}
		}
		else //if the cell will get infected through in this time step
		{
			prev = c;
		}
	}

	//###THIS GETS PARTICULARLY PROBLEMATIC FOR STRONG DIFFERENCES IN PARTICLE NUMBERS PER NODE
	//###AND NODES DO NOT ALWAYS HAVE THE SAME SIZE!

	cells.data.totalRXfrac = totalRXsum / ((Real)cells.data.zFac * (Real)cells.data.zPer);
	gen.delta = MIN(gen.ASPIREDFILLPERSTEP * cells.data.size, .2*cells.data.size/V);
	// cout << Tcells->shortesttNucleationtime << endl; for debug
		return somethingHasMoved;
}

Real caHdl::currentCellTemperature(cellP c, Real time, long timeKey)  //This can be somewhere else. A class for temperature management is a better option
{
	QUICKASSERT(cellExists(c));

	short ix = c->ix;
	short iy = c->iy;
	short iz = c->iz;


	long prevKey = 0;
	long nextKey = Tcells->T_timer.size - 1;

	TCellP Tcell = TCellOfActiveCell(c);

	if (timeKey < Tcells->T_timer.size - 1) nextKey = timeKey + 1;



	Real Ti = Tcell->interpolateTemperatureSpace(ix, iy, iz, time, timeKey, nextKey);

	return Ti;
}

void caHdl::updateFullRX(void)
{
	cellP c, nextCell;

	//currently //compromise //originalfactors //seemed optimized factors but are not //original factors for velocity fgeo 1 in <100>, 1/2^0.5 <110>, 1/3^0.5 <111>
	const Real overFac1 = 0.0682;// 1;// 0.0682(orig); //1;//  0.0682; //0.04946756; //0.03738274; //0.0523448297609848;  //###MK20121211, 26NN, 8 face, 12 edge and 6 diag neighbors each one getting a portion of the overshoot, all other overshoot is discarded
	const Real overFac2 =  0.0542;// 0.0542 (orig) //0.7071;; //0.0387; //0.03949696; //0.04059546; //0.0370133840840478;
	const Real overFac3 =   0.0484;// 0.0484;(orig)//0.577; //0.0307; //0.02865390; //0.03606975; //0.0302213015531897;


	for (c = cells.first[0]; c; c = nextCell)
	{

		nextCell = c->next;

		Real overshoot = c->rxFrac - 1.0;            //Compute overshoot when cell is full Rx //Berechne Overshoot, wenn Zelle voll Rx
		Real overshoot1 = overshoot * overFac1;
		Real overshoot2 = overshoot * overFac2;
		Real overshoot3 = overshoot * overFac3;

		//check: instead of defGrain Tcell
		//check: set time for nucleation for infected cell

		grainP rxg = grainOfCell(c);		//rx-grain
		//grainP dfg = defGrainOfCell( c );	//deformed grain 

		//TCellP TCell = TCellOfActiveCell(c);
		//long indexTCellInfector = TCell->id;

		int indexInfector = rxg->oriIndex;			//int indexInfector = infectorCell->grain->oriIndex;
		//int indexInfected = dfg->oriIndex;			//int indexInfected = c->grain->oriIndex;
		int indexInfected = 0;					//only the first grain is the melt

		grains.orientations[indexInfector]->count++;  //cell-based discrete texture update shortly before cell retirement
		grains.orientations[indexInfected]->count--;

		rxg->CellCount++;                   //here: Count up the number of 100 % Rx cells of the grain by 1 //hier:Anzahl der 100%Rx-Cells des Korns um 1 hochzaehlen

		//NOW c becomes itself an infector for new victims, we view the neighborhood as the infector sees it but assign the values in the view of the victim
		// for AM we have to check of it is write time to infect it or not.

		infectCell(c, -1, 0, 0, 13, overshoot1); //25 is opposite to 0, each infector char is so to say inverted at the centre 0,0,0
		infectCell(c, 0, -1, 0, 15, overshoot1);  //those dx,dy,dz coordinates are directly passed through to sendtointerface
		infectCell(c, 0, 0, -1, 21, overshoot1);
		infectCell(c, 0, 0, +1, 4, overshoot1);
		infectCell(c, 0, +1, 0, 10, overshoot1);
		infectCell(c, +1, 0, 0, 12, overshoot1);

		infectCell(c, +1, +1, 0, 9, overshoot2);	//xy (edge)
		infectCell(c, -1, +1, 0, 11, overshoot2);
		infectCell(c, +1, -1, 0, 14, overshoot2);
		infectCell(c, -1, -1, 0, 16, overshoot2);

		infectCell(c, 0, +1, +1, 1, overshoot2);	//yz(edge)
		infectCell(c, 0, -1, +1, 7, overshoot2);
		infectCell(c, 0, +1, -1, 18, overshoot2);
		infectCell(c, 0, -1, -1, 24, overshoot2);

		infectCell(c, +1, 0, +1, 3, overshoot2);	//xz(edge)
		infectCell(c, -1, 0, +1, 5, overshoot2);
		infectCell(c, +1, 0, -1, 20, overshoot2);
		infectCell(c, -1, 0, -1, 22, overshoot2);

		infectCell(c, +1, +1, +1, 0, overshoot3);	//xyz(diag)
		infectCell(c, -1, +1, +1, 2, overshoot3);
		infectCell(c, +1, -1, +1, 6, overshoot3);
		infectCell(c, -1, -1, +1, 8, overshoot3);
		infectCell(c, +1, +1, -1, 17, overshoot3);
		infectCell(c, -1, +1, -1, 19, overshoot3);
		infectCell(c, +1, -1, -1, 23, overshoot3);
		infectCell(c, -1, -1, -1, 25, overshoot3);

		retireCell(c);                      //cell for recycling
	}
	long i;
	for (i = 0; i < TIMESLOTCOUNT; i++)			cells.first[i] = cells.first[i + 1];
	cells.first[TIMESLOTCOUNT] = NULL;
}

cellP caHdl::infectCell(cellP infector, long dx, long dy, long dz, int infectorvictimindex, Real overshoot)  //explicit cast of int to char 0 <= infectorvictimindex is only possible because <= 25 
{

	long iz = infector->iz + dz;  //many bytes transferred in the header
	long iy = infector->iy + dy;  //this is the global position of the cell which is to be infected
	long ix = infector->ix + dx;

	if (handlePeriodicity(&ix, &iy, &iz) == NOTHINGTOINFECT) return NULL;

	if (COMM_ENABLED & protocol)
	{
		int nextNode = parallelHdl->cellNode(ix, iy, iz);
		int mynode = myNode();

		if (nextNode != mynode)
		{
			sendToInterface(infector, dx, dy, dz, nextNode, overshoot); //###LB: overshoot not yet used
			return NULL;
		}
	}

	//SO IF THE INFECTION PROCESS IS TO BE PERFORMED IN LOCAL NODE

	cellP* ncp = &cells.all[(iz - z0) * cells.data.zFac + (iy - y0) * cells.data.xPer + (ix - x0)];
	cellP nc = *ncp;

	//if( cellExists(nc) )	return nc; //###LB: Changed

	if (cellExists(nc))
	{
		if (!InsufficientMelting(nc) && !isParticle(nc) && isCellRecrystallized(nc) && (grainOfCell(nc) != grainOfCell(infector)))
		{
			cell cb;
			cb.ix = (short)ix;
			cb.iy = (short)iy;
			cb.iz = (short)iz;
			cb.dataHolder = nc->dataHolder;

			/*long s = SQR( dx ) + SQR( dy ) + SQR( dz );

			if( s == 1 )
				new cellsBoundary( &boundaries, &cb, infector );*/ //Can be used to track boundaries
		}
		return NULL;
	}

	unsigned char rhoFac_RFU = ((long)nc & 0x7FF) >> 3; //Mask with 0x7FF is not really necessary because the value is cast in a char but...

	long defGrainPos = grainArrayPosition(nc);		//It returns the index of the melt (0); it is correct
	long rxGrainPos = posOfGrainOfCell(infector);	//rx grain

	Real Temperature = cellTemperature(nc) * TempFactor; //Check: This is not needed anymore

	TCellP TCellc = getTCellofCell(nc);
	long TCellID = TCellc->id;

	//Real nucTemperature = TCellc->determineNucleationTime(ix, iy, iz);
	Real nucTemperature = Tcells->getCellNucleationTime(nc);

	nc = getCellMemory();

	*ncp = nc;

	nc->ix = (short)ix;
	nc->iy = (short)iy;
	nc->iz = (short)iz;

	grainP rxg = grains.all[rxGrainPos];

	if (ix > rxg->mySpace.xright)	rxg->mySpace.xright = (short)ix;
	if (ix < rxg->mySpace.xleft)	rxg->mySpace.xleft = (short)ix;

	if (iy > rxg->mySpace.yright)	rxg->mySpace.yright = (short)iy;
	if (iy < rxg->mySpace.yleft)	rxg->mySpace.yleft = (short)iy;

	if (iz > rxg->mySpace.zright)	rxg->mySpace.zright = (short)iz;
	if (iz < rxg->mySpace.zleft)	rxg->mySpace.zleft = (short)iz;

	nc->rxFrac = overshoot;

	nc->dataHolder = 0;

	nc->rhoFac_RFU = rhoFac_RFU;

	char infidx = (char)infectorvictimindex;
	//uint_fast64_t dgidx = (uint_fast64_t) defGrainPos << 43;
	uint_fast64_t TCidx = (uint_fast64_t)TCellID << 43;
	uint_fast64_t rxgidx = (uint_fast64_t)rxGrainPos << 22;

	nc->dataHolder |= infidx;
	nc->dataHolder |= rxgidx;

	//nc->dataHolder |= dgidx;

	nc->next = cells.first[1];
	cells.first[1] = nc;

	uint_fast64_t timeDiscrete = 0xFFFF;

	if((long)nucTemperature != LIQUID && (long)nucTemperature != NO_MELT) //as:changed
		timeDiscrete = (uint_fast64_t)(nucTemperature * Tcells->timeEncodeFactor);

	timeDiscrete <<= 5;

	nc->dataHolder |= timeDiscrete;
	nc->dataHolder |= TCidx;

	nc->T = Temperature;

	//grainP dg = defGrainOfCell(nc); //deformed grain of cell

	//dg->CellCount--;

	cells.data.count++; //cells currently active in main loop, even if camouflaged for jump timesteps by ASPIREDFILLPERSTEP
	return nc;
}

int caHdl::sendToInterface(cellP infector, long dx, long dy, long dz, int nextNode, Real overshoot)  //those dx,dy,dz are the relative coordinates placing the automaton coordinate system into the infector
{
	int ix = infector->ix; //global position of the cell which infects the new cell
	int iy = infector->iy;
	int iz = infector->iz;

	int dxyzdxyz = dx * dx + dy * dy + dz * dz;
	QUICKASSERT(dxyzdxyz == 1 || dxyzdxyz == 2 || dxyzdxyz == 3); //only for debugging purposes

	if (dxyzdxyz == 1) { //rejects all 110 and 111 infection attempts

		if (dz == -1) {
			QUICKASSERT(nextNode == nodeFace[C_FRONT]->nextNode);
			nodeFace[C_FRONT]->addCell100(infector, ix, iy);
			somethingToSend |= nodeFace[C_FRONT]->somethingToSend;
			return 1;
		}

		if (dz == 1) {
			QUICKASSERT(nextNode == nodeFace[C_REAR]->nextNode);
			nodeFace[C_REAR]->addCell100(infector, ix, iy);
			somethingToSend |= nodeFace[C_REAR]->somethingToSend;
			return 1;
		}

		if (dy == -1) {
			QUICKASSERT(nextNode == nodeFace[C_BOTTOM]->nextNode);
			nodeFace[C_BOTTOM]->addCell100(infector, ix, iz);
			somethingToSend |= nodeFace[C_BOTTOM]->somethingToSend;
			return 1;
		}

		if (dy == 1) {
			QUICKASSERT(nextNode == nodeFace[C_TOP]->nextNode);
			nodeFace[C_TOP]->addCell100(infector, ix, iz);
			somethingToSend |= nodeFace[C_TOP]->somethingToSend;
			return 1;
		}

		if (dx == -1) {
			QUICKASSERT(nextNode == nodeFace[C_LEFT]->nextNode);
			nodeFace[C_LEFT]->addCell100(infector, iz, iy);
			somethingToSend |= nodeFace[C_LEFT]->somethingToSend;
			return 1;
		}

		if (dx == 1) {
			QUICKASSERT(nextNode == nodeFace[C_RIGHT]->nextNode);
			nodeFace[C_RIGHT]->addCell100(infector, iz, iy);
			somethingToSend |= nodeFace[C_RIGHT]->somethingToSend;
			return 1;
		}

	}
	return 0;
}

char caHdl::handlePeriodicity(long* x, long* y, long* z)
{
	long ix = *x;
	long iy = *y;
	long iz = *z;

	long xglobal = box.xPer;
	long yglobal = box.yPer;
	long zglobal = box.zPer;

	switch (boxPer & nodePer)
	{
	case (BOXPER & NON_NODE_PER):

		if (iz < 0)		iz += zglobal;
		if (iz >= zglobal)	iz -= zglobal;

		if (iy < 0)		iy += yglobal;
		if (iy >= yglobal)	iy -= yglobal;

		if (ix < 0)		ix += xglobal;
		if (ix >= xglobal)	ix -= xglobal;
		break;

	case (NONBOXPER & NODE_PER):

		if (iz < z0)				iz += zPer;
		if (iz >= (z0 + zPer))	iz -= zPer;

		if (iy < y0)				iy += yPer;
		if (iy >= (y0 + yPer))	iy -= yPer;

		if (ix < x0)				ix += xPer;
		if (ix >= (x0 + xPer))	ix -= xPer;
		break;

	case(NONBOXPER & NON_NODE_PER):

		if (iz < 0)					return FAILURE;
		if (iz >= zglobal)				return FAILURE;
		if (iy < 0)					return FAILURE;
		if (iy >= yglobal)				return FAILURE;
		if (ix < 0)					return FAILURE;
		if (ix >= xglobal)				return FAILURE;

		if (protocol == NOCOMM)
		{
			if (iz < z0)				return FAILURE;
			if (iz >= (z0 + zPer))	return FAILURE;
			if (iy < y0)				return FAILURE;
			if (iy >= (y0 + yPer))	return FAILURE;
			if (ix < x0)				return FAILURE;
			if (ix >= (x0 + xPer))	return FAILURE;
		}
		break;
	}

	*x = ix;
	*y = iy;
	*z = iz;

	return SUCCESS;
}

cellP caHdl::getCellMemory(void)
{
	cellP nc = cells.firstRecycled;

	if (!nc)
	{
		long i;
#define BLOCKSIZE (32768-32)/sizeof(cell) //####not portable!

		nc = (cellP)malloc(sizeof(cell) * BLOCKSIZE);
		if (!nc) err->reportError(ERRTXT("Cannot allocate more memory"));

		for (i = 0; i < BLOCKSIZE; i++)
		{
			cellP c = &nc[i];
			c->next = &nc[i + 1];
		}
		nc[BLOCKSIZE - 1].next = NULL; //Concatenate all cells in this block
		//as malloc allocates a continous block on the virtual address space but not necessarily continous in physical address space
		//this procedure limits page misses: EXACTLY
		cells.firstRecycled = nc;
		cells.countRecycled += BLOCKSIZE;
	}

	//###utilize this construct by initializing the next cell in this VRM continuous chunk as firstRecycled?
	cells.firstRecycled = nc->next;
	cells.countRecycled--;
	return nc;
}

void caHdl::retireCell(cellP c)
{

	cells.data.count--;   //The number of cells is reduced by one
	cells.data.retiredCount++;

	grainP gr = grainOfCell(c); //Recrystallized cell belongs to a recrystallized grain

	//grainP gr = c->grain;    //a cell must become a zombie and receives grain pointers from the cell//eine Zelle muss zum Zombie werden und erhaelt grain zeiger von der Zelle

	if (!gr->zombie)    //if there is no zombie yet; Zombie contains belonging to the grain
	{
		c->rxFrac = 1000.0;   //rx fract of the notched cell becomes 1//rx fract der ausgeklinkten Zelle wird 1
		gr->zombie = c;    //c becomes a zombie
		return;
	}

	*getCellPerP(c->ix, c->iy, c->iz) = gr->zombie;
	freeCellMemory(c);
}

void caHdl::freeCellMemory(cellP c)
{//Cell memory is reused later
	c->next = cells.firstRecycled;
	cells.firstRecycled = c;
	cells.countRecycled++;
}

void caHdl::sendReceive(void)
{
	int currentNode = myNode();

	for (int sender = 0; sender < totalNodes; sender++)
	{
		if (sender == currentNode)
		{
			world::sendReceive();
			somethingToSend = NOTHING;
		}
		else
		{
			receiveSend(sender);
		}
	}

	/*	parallelDataP nodeData = parallelHdl->all[currentNode];

		ostringstream message;

		int ixn = nodeData->ixn;
		int iyn = nodeData->iyn;
		int izn = nodeData->izn;

		//parity
		bool xpar = EVEN;
		bool ypar = EVEN;
		bool zpar = EVEN;

		if( ixn % 2 ) xpar = ODD;
		if( iyn % 2 ) ypar = ODD;
		if( izn % 2 ) zpar = ODD;

		bool nodeState = !XOR(XOR(xpar,ypar),zpar);

		int st=0;
		MPI_Status status;
		int stSend = (int) nodeState;

		message << " Node: " << myNode() << " ns1: " << nodeState << endl;

		for( int i=1;i<totalNodes;i++ )
		{
			if( myNode()==0 )
			{
				MPI_Recv( &st,1,MPI_INT,i,0,MPI_COMM_WORLD, &status );
				message << " Node: " << i << " ns1: " << st << endl;
			}
			else
				MPI_Send( &stSend,1,MPI_INT,0,0,MPI_COMM_WORLD );
		}

		write( message.str(), &logFile );

		if( nodeState )
		{
			interfaceHdl::sendReceive();
			somethingToSend = NOTHING;
		}else
		{
			receiveSend( sender );
		}*/

}

void caHdl::receiveSend(int sender)
{

	parallelDataP sendingNode = this->parallelHdl->all[sender];
	//###20121228 myNode()
	if (sendingNode->neighbourNodes[C_RIGHT] == myNode())	world::receiveSend(sender);
	if (sendingNode->neighbourNodes[C_FRONT] == myNode())	world::receiveSend(sender);
	if (sendingNode->neighbourNodes[C_TOP] == myNode())	world::receiveSend(sender);
	if (sendingNode->neighbourNodes[C_BOTTOM] == myNode())	world::receiveSend(sender);
	if (sendingNode->neighbourNodes[C_REAR] == myNode())	world::receiveSend(sender);
	if (sendingNode->neighbourNodes[C_LEFT] == myNode())	world::receiveSend(sender);

	/*
	//###20121228done <110> and <111> sender
	if( sendingNode->neighbourNodes[C_FRONTTOP]			== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_FRONTBOTTOM]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_FRONTLEFT]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_FRONTRIGHT]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARTOP]			== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARBOTTOM]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARLEFT]			== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARRIGHT]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_RIGHTTOP]			== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_RIGHTBOTTOM]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_LEFTTOP]			== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_LEFTBOTTOM]		== myNode() )	world::receiveSend( sender );

	if( sendingNode->neighbourNodes[C_FRONTTOPLEFT]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_FRONTTOPRIGHT]	== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_FRONTBOTTOMLEFT]	== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_FRONTBOTTOMRIGHT]	== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARTOPLEFT]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARTOPRIGHT]		== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARBOTTOMLEFT]	== myNode() )	world::receiveSend( sender );
	if( sendingNode->neighbourNodes[C_REARBOTTOMRIGHT]	== myNode() )	world::receiveSend( sender );
	*/
}

void caHdl::receiveSend(void)
{
	//###20121228done
	for (int direction = 0; direction < 6; direction++)
		nodeFace[direction]->sendDataRecvData();

	world::receiveSend();
}

void caHdl::distributeReceivedCells(void)
{
	for (int i = 0; i < SIXFACES; i++) {
		if (nodeFace[i]->somethingReceived) {
			if (protocol == DETERMINISTIC) nodeFace[i]->distributeImportedCells100();
			//###20121228 statistical mode not implemented yet
			if (protocol == STATISTICAL)	nodeFace[i]->distributeImportedCells(insertionSite);
		}
	}

	somethingReceived = NOTHING;
}

cellP caHdl::addTransferredCell(grainP g, long ix, long iy, long iz) // Used by interface.cpp
{
	if (g->DeformGrainIndex != -1) err->reportError(ERRTXT("It must be a non-deformed grain"));

	cell pseudoNb;                     //construct needed to infect cell, that is the germ cell in the recrystallizing grain// benoetigtes Konstrukt um Zelle zu infizieren, das ist die Keimzelle im rekristallisierenden Korn
	pseudoNb.ix = (short)ix;          // what happens if ix is ​​outside the short limits ### portability? When would you need 32767 ^ 3 cells?//was passiert wenn ix außerhalb der short Grenzen ist ###portability? When would you need 32767^3 cells?
	pseudoNb.iy = (short)iy;			//what happens if I am outside the short limits ### portability? When would you need 32767^3 cells?
	pseudoNb.iz = (short)iz;

	uint_fast64_t gpos = (uint_fast64_t)g->arrayPos;

	pseudoNb.dataHolder = (0 | (gpos << 22));

	cellP nucCell = infectCell(&pseudoNb, 0, 0, 0, 26, 0.0);

	if (nucCell == NOTNEWCELLCREATED) return NULL;

	QUICKASSERT(nucCell);

	return nucCell;
}

void caHdl::endParallelProcesses(void)
{
	MPI_Finalize();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// These were used for writing binary file in first version of microstructure generator////////////
////// Not needed for Additive Manufacturing simulations.//////////////////////

void caHdl::addOrientation(Real phi1, Real PHI, Real phi2) // No use.
{
	grainP nucGrain = new grain(&(this->grains), phi1, PHI, phi2, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 1, -1);
	QUICKASSERT(nucGrain);
}

void caHdl::addNucleus(cellsToNucleate* cn, short ix, short iy, short iz, grainP rxg)  //not callled
{
	nucleiP nuc = (nucleiP)calloc(1, sizeof(nuclei));

	if (!nuc) err->reportError(ERRTXT("Cannot allocate more memory"));

	nuc->grainIdx = rxg->arrayPos;
	nuc->ix = ix;
	nuc->iy = iy;
	nuc->iz = iz;

	if (cn->first) cn->first->prev = nuc;
	nuc->next = cn->first;
	nuc->prev = NULL;
	cn->first = nuc;
	if (!(cn->last)) cn->last = nuc;

	cn->count++;
	cn->nucleateFlag = 0;

	rxg->mySpace.xleft = rxg->mySpace.xright = ix;
	rxg->mySpace.yleft = rxg->mySpace.yright = iy;
	rxg->mySpace.zleft = rxg->mySpace.zright = iz;

	rxg->CellCount = 0;
	rxg->initialCellNumber = 0;
}

void caHdl::saveMicrostructure(const char* filename) //not called
{
	ostringstream message;

	MPI_Offset offset = 0;

	char* fname = new char[BUFSIZ];

	strcpy(fname, filename);

	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &theFile);

	write("Writing file header", &logFile);
	writeHeader(&offset);
	write("Writing orientations", &logFile);
	writeOrientations(&offset);
	write("Writing grains", &logFile);
	nucleiPacket** packets = writeGrains(&offset);
	write("Writing nuclei", &logFile);
	writeNuclei(&offset, &packets);
	write("Writing cells", &logFile);
	writeCells(&offset);
	write("Closing file", &logFile);
	MPI_File_close(&theFile);
	theFile = NULL;

	delete fname;
}

void caHdl::writeHeader(MPI_Offset* offset)		//###LB: Last modified: 22.07.2013
{
	MPI_Offset localOffset = *offset;

	int rank = myNode();

	long identifier[4] = { IDTAG0, IDTAG1, IDTAG2, IDTAG3 };
	long dimensions[4] = { (long)totalNodes, box.xPer, box.yPer, box.zPer };
	char periodicity[4] = { nodePer, protocol, commType, insertionSite };
	double boxProperties[6] = { box.Volume, box.size, cells.data.size, xsize, ysize, zsize };  //###LB: box properties added to output

	if (rank == MASTER)
	{
		MPI_File_write_at(theFile, localOffset, identifier, 4, MPI_LONG, MPI_STATUS_IGNORE);

		localOffset += 4 * sizeof(long);

		MPI_File_write_at(theFile, localOffset, dimensions, 4, MPI_LONG, MPI_STATUS_IGNORE);

		localOffset += 4 * sizeof(long);

		MPI_File_write_at(theFile, localOffset, boxProperties, 6, MPI_DOUBLE, MPI_STATUS_IGNORE);

		localOffset += 6 * sizeof(double);

		MPI_File_write_at(theFile, localOffset, periodicity, 4, MPI_CHAR, MPI_STATUS_IGNORE);

		localOffset += 4 * sizeof(char);
	}
	else
	{
		localOffset += (8 * sizeof(long) + 4 * sizeof(char) + 6 * sizeof(double));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	localOffset += 6 * sizeof(long) * rank;

	long localDim[6] = { x0, y0, z0, xPer, yPer, zPer };

	MPI_File_write_at(theFile, localOffset, localDim, 6, MPI_LONG, MPI_STATUS_IGNORE);

	localOffset += 6 * sizeof(long) * (totalNodes - rank);

	*offset = localOffset;

}

void caHdl::writeOrientations(MPI_Offset* offset)
{
	MPI_Offset localOffset = *offset;
	if (rank == MASTER)
	{
		long Nori = grains.noris;

		MPI_Offset tof = localOffset;

		MPI_File_write_at(theFile, localOffset, &Nori, 1, MPI_LONG, MPI_STATUS_IGNORE);

		localOffset += sizeof(long);

		exportOri* oris = (exportOri*)calloc(Nori, sizeof(exportOri));
		if (!oris) err->reportError(ERRTXT("Cannot allocate more memory"));

		for (long i = 0; i < Nori; i++)
		{
			oris[i].count = grains.orientations[i]->count;
			oris[i].PHI = grains.orientations[i]->PHI;
			oris[i].phi1 = grains.orientations[i]->phi1;
			oris[i].phi2 = grains.orientations[i]->phi2;
			oris[i].q0 = grains.orientations[i]->q0;
			oris[i].q1 = grains.orientations[i]->q1;
			oris[i].q2 = grains.orientations[i]->q2;
			oris[i].q3 = grains.orientations[i]->q3;
		}
		MPI_File_write_at(theFile, localOffset, oris, Nori, MPI_IO_Ori_Type, MPI_STATUS_IGNORE);

		localOffset += Nori * sizeof(exportOri);

		free(oris);
	}
	MPI_Bcast(&localOffset, 1, MPI_LONG_LONG, MASTER, MPI_COMM_WORLD);

	*offset = localOffset;
}

nucleiPacket** caHdl::writeGrains(MPI_Offset* displ)
{
	MPI_Offset localOffset = *displ;

	MPI_IO_GrainInfo* grainData = NULL;
	nucleiPacket** packets = NULL;

	if (myNode() == MASTER)
	{
		long arraysSizes[2] = { grains.Nuc, grains.NucMG };
		long size = grains.Nuc + grains.NucMG;

		ostringstream message;

		message << "NumberOfNuclei:" << grains.NucMG << " NumberOfDefGrains:" << grains.Nuc << endl;

		write(message.str(), &logFile);

		grainData = (MPI_IO_GrainInfo*)calloc(size, sizeof(MPI_IO_GrainInfo));
		if (!grainData) err->reportError(ERRTXT("Cannot allocate more memory"));

		packets = (nucleiPacket**)calloc(4 * grains.Nuc, sizeof(nucleiPacket)); //###LB: deformed grains
		if (!packets) err->reportError(ERRTXT("Cannot allocate more memory"));

		for (long i = 0; i < size; i++)
		{
			grainP g = grains.all[i + grains.defGrainsCount];
			prepareGrainForBroadcast(i + grains.defGrainsCount, &grainData[i]);
			grainData[i].arrayPos -= grains.defGrainsCount;

			if ((i + grains.defGrainsCount) >= (grains.defGrainsCount + grains.Nuc)) //###LB: Because nuclei grains do not have nuclei.
				continue;

			if (g->GBnuclei.count == 0) packets[4 * i] = NULL;
			else
			{
				nucleiPacket* gbpack = (nucleiPacket*)malloc(sizeof(nucleiPacket));
				if (!gbpack) err->reportError(ERRTXT("Cannot allocate more memory"));

				gbpack->count = g->GBnuclei.count;

				gbpack->nucList = (exportNucleus**)calloc(gbpack->count, sizeof(exportNucleus*));
				if (!gbpack->nucList) err->reportError(ERRTXT("Cannot allocate more memory"));

				int counter = 0;

				nucleiP ncnext = NULL;

				for (nucleiP ncp = g->GBnuclei.first; ncp; ncp = ncnext)
				{
					ncnext = ncp->next;

					exportNucleus* nuc = (exportNucleus*)calloc(1, sizeof(exportNucleus));
					if (!nuc) err->reportError(ERRTXT("Cannot allocate more memory"));

					nuc->x = ncp->ix;
					nuc->y = ncp->iy;
					nuc->z = ncp->iz;
					gbpack->nucList[counter] = nuc;
					counter++;
					free(ncp);
				}

				packets[4 * i] = gbpack;

			}

			if (g->SBnuclei.count == 0) packets[4 * i + 1] = NULL;
			else
			{
				nucleiPacket* sbpack = (nucleiPacket*)malloc(sizeof(nucleiPacket));
				if (!sbpack) err->reportError(ERRTXT("Cannot allocate more memory"));

				sbpack->count = g->SBnuclei.count;

				sbpack->nucList = (exportNucleus**)calloc(sbpack->count, sizeof(exportNucleus*));
				if (!sbpack->nucList) err->reportError(ERRTXT("Cannot allocate more memory"));

				int counter = 0;

				nucleiP ncnext = NULL;

				for (nucleiP ncp = g->SBnuclei.first; ncp; ncp = ncnext)
				{
					ncnext = ncp->next;

					exportNucleus* nuc = (exportNucleus*)calloc(1, sizeof(exportNucleus));
					if (!nuc) err->reportError(ERRTXT("Cannot allocate more memory"));

					nuc->ig = ncp->grainIdx - grains.defGrainsCount;
					nuc->x = ncp->ix;
					nuc->y = ncp->iy;
					nuc->z = ncp->iz;
					sbpack->nucList[counter] = nuc;
					counter++;
					free(ncp);
				}

				packets[4 * i + 1] = sbpack;
			}

			if (g->TBnuclei.count == 0) packets[4 * i + 2] = NULL;
			else
			{
				nucleiPacket* tbpack = (nucleiPacket*)malloc(sizeof(nucleiPacket));
				if (!tbpack) err->reportError(ERRTXT("Cannot allocate more memory"));

				tbpack->count = g->TBnuclei.count;

				tbpack->nucList = (exportNucleus**)calloc(tbpack->count, sizeof(exportNucleus*));
				if (!tbpack->nucList) err->reportError(ERRTXT("Cannot allocate more memory"));

				int counter = 0;

				nucleiP ncnext = NULL;

				for (nucleiP ncp = g->TBnuclei.first; ncp; ncp = ncnext)
				{
					ncnext = ncp->next;

					exportNucleus* nuc = (exportNucleus*)calloc(1, sizeof(exportNucleus));
					if (!nuc) err->reportError(ERRTXT("Cannot allocate more memory"));

					nuc->ig = ncp->grainIdx - grains.defGrainsCount;
					nuc->x = ncp->ix;
					nuc->y = ncp->iy;
					nuc->z = ncp->iz;
					tbpack->nucList[counter] = nuc;
					counter++;
					free(ncp);
				}

				packets[4 * i + 2] = tbpack;
			}

			if (g->PSNnuclei.count == 0) packets[4 * i + 3] = NULL;
			else
			{
				nucleiPacket* psnpack = (nucleiPacket*)malloc(sizeof(nucleiPacket));
				if (!psnpack) err->reportError(ERRTXT("Cannot allocate more memory"));

				psnpack->count = g->PSNnuclei.count;

				psnpack->nucList = (exportNucleus**)calloc(psnpack->count, sizeof(exportNucleus*));
				if (!psnpack) err->reportError(ERRTXT("Cannot allocate more memory"));

				int counter = 0;

				nucleiP ncnext = NULL;

				for (nucleiP ncp = g->PSNnuclei.first; ncp; ncp = ncnext)
				{
					ncnext = ncp->next;

					exportNucleus* nuc = (exportNucleus*)calloc(1, sizeof(exportNucleus));
					if (!nuc) err->reportError(ERRTXT("Cannot allocate more memory"));

					nuc->ig = ncp->grainIdx - grains.defGrainsCount;
					nuc->x = ncp->ix;
					nuc->y = ncp->iy;
					nuc->z = ncp->iz;
					psnpack->nucList[counter] = nuc;
					counter++;
					free(ncp);
				}

				packets[4 * i + 3] = psnpack;
			}
		}

		MPI_Offset gsize = sizeof(MPI_IO_GrainInfo);

		MPI_File_write_at(theFile, localOffset, &arraysSizes, 2, MPI_LONG, MPI_STATUS_IGNORE);

		localOffset += 2 * sizeof(long);

		MPI_File_write_at(theFile, localOffset, grainData, size, MPI_IO_GrainInfo_Type, MPI_STATUS_IGNORE);

		localOffset += size * gsize;

		free(grainData);
	}

	for (long i = 0; i < grains.noris; i++)
	{
		free(grains.orientations[i]);
		grains.orientations[i] = NULL;
	}

	MPI_Bcast(&localOffset, 1, MPI_LONG_LONG, MASTER, MPI_COMM_WORLD);
	*displ = localOffset;
	return packets;
}

void caHdl::writeNuclei(MPI_Offset* disp, nucleiPacket*** packP)
{
	MPI_Offset localOffset = *disp;

	nucleiPacket** packets = *packP;

	if (rank == MASTER)
	{
		for (long i = 0; i < 4 * grains.Nuc; i += 4)
		{
			for (int j = 0; j < 4; j++)
			{

				nucleiPacket* p = packets[i + j];
				long arraySize = 0;

				if (p) arraySize = p->count;

				MPI_File_write_at(theFile, localOffset, &arraySize, 1, MPI_LONG, MPI_STATUS_IGNORE);
				localOffset += sizeof(long);

				long counter = 0;

				while (counter < arraySize)
				{
					exportNucleus* toSaveP = packets[i + j]->nucList[counter];

					MPI_File_write_at(theFile, localOffset, toSaveP, 1, MPI_IO_Nucleus_Type, MPI_STATUS_IGNORE);

					localOffset += sizeof(exportNucleus);

					free(toSaveP);

					packets[i + j]->nucList[counter] = NULL;

					counter++;
				}
				if (packets[i + j])
				{
					free(packets[i + j]->nucList);
					free(packets[i + j]);
					packets[i + j] = NULL;
				}
			}
		}
	}
	if (packets)	free(packets);
	MPI_Bcast(&localOffset, 1, MPI_LONG_LONG, MASTER, MPI_COMM_WORLD);
	*disp = localOffset;
}

/*void caHdl::writeParticles(MPI_Offset* disp)
{
	MPI_Offset localOffset = *disp;
	long pcount = 0;

	if( rank == MASTER )
	{

		if( particles )
			pcount = particles->count;

		MPI_File_write_at( theFile, localOffset, &pcount, 1, MPI_LONG, MPI_STATUS_IGNORE );
		localOffset+= sizeof( long );

		if( !particles )
		{
			*disp = localOffset;
			MPI_Bcast( &localOffset, 1, MPI_LONG_LONG, MASTER, MPI_COMM_WORLD );
			return;
		}

		Real dxyz[3] = { particles->ZONEDX, particles->ZONEDY, particles->ZONEDZ };

		MPI_File_write_at( theFile, localOffset, dxyz, 3, MPI_DOUBLE, MPI_STATUS_IGNORE );
		localOffset += 3 * sizeof(Real);

		exportParticle * pa = (exportParticle *) calloc( pcount, sizeof(exportParticle) );
		if( !pa ) err->reportError(ERRTXT("Cannot allocate more memory"));

		for( int i=0; i < pcount; i++ )
		{
			particleP p = particles->all[i];

			pa[i].count = p->cellCount;
			pa[i].pos = p->arrayPos;

			pa[i].xl = p->leftX;	pa[i].yl = p->leftY;	pa[i].zl = p->leftZ;
			pa[i].xr = p->rightX;	pa[i].yr = p->rightY;	pa[i].zr = p->rightZ;
		}

		MPI_File_write_at( theFile, localOffset, pa, pcount, MPI_IO_Particle_Type, MPI_STATUS_IGNORE );

		localOffset += pcount * sizeof(exportParticle);
		free( pa );
	}
	MPI_Bcast( &localOffset, 1, MPI_LONG_LONG, MASTER, MPI_COMM_WORLD );
	*disp = localOffset;
}*/

void caHdl::writeCells(MPI_Offset* disp)
{
	long size;
	MPI_Offset localOffset = *disp;

	long* compCells = compressNodeInformation(&size);

	long* compSize = NULL;

	size *= 2;

	if (rank == MASTER) compSize = (long*)calloc(totalNodes, sizeof(long));

	MPI_Gather(&size, 1, MPI_LONG, compSize, 1, MPI_LONG, MASTER, MPI_COMM_WORLD);

	if (rank == MASTER) MPI_File_write_at(theFile, localOffset, compSize, totalNodes, MPI_LONG, MPI_STATUS_IGNORE);

	localOffset += totalNodes * sizeof(long);

	long long sizeBcast = size;

	long long longOffset = 0;

	int rankCount = 0;

	while (rankCount < (totalNodes - 1))
	{
		MPI_Bcast(&sizeBcast, 1, MPI_LONG_LONG, rankCount, MPI_COMM_WORLD);

		if (rank > rankCount) longOffset += sizeBcast;

		sizeBcast = size;

		rankCount++;
	}

	localOffset += sizeof(long) * longOffset;	//For master, longOffset = 0;

	MPI_File_write_at(theFile, localOffset, compCells, size, MPI_LONG, MPI_STATUS_IGNORE);

	free(compCells);
}

long* caHdl::compressNodeInformation(long* size) //###LB: RLE
{
	size_t arrayLength = 5003 * 2;//###LB: Twice a prime

	long* temp = (long*)calloc(arrayLength, sizeof(long));
	if (!temp) err->reportError(ERRTXT("Cannot allocate more memory"));

	long count = 0;

	cellP c = cells.all[0];
	long c0 = (long)c;

	temp[0] = c0;
	temp[1] = 0;

	for (long i = 0; i < xPer * yPer * zPer; i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);

		cellP cl = cells.all[i];

		long ccurrent = (long)cl;

		if (ccurrent == c0)
		{
			temp[count + 1]++;
		}
		else
		{
			if (((count + 2) / 2) % (arrayLength / 2) == 0)
			{
				size_t newLength = ((count + (long)2) / arrayLength + 1) * arrayLength;

				long* newArray = (long*)calloc(newLength, sizeof(long));

				if (!newArray) err->reportError(ERRTXT("Cannot allocate more memory"));

				long* toDel = temp;

				memcpy(newArray, temp, (count + (long)2) * sizeof(long));

				free(toDel);
				temp = newArray;
			}

			if ((temp[count] & 0x7) != 0x2)  //###LB: If cell is not a particle the grain array position will be displaced by the defgrainscount offset
			{
				long ct = temp[count];
				long gi = ct >> 11;
				gi -= grains.defGrainsCount;

				if (gi < 0)
					err->reportError(ERRTXT("The grain index cannot be lower than 0"));

				gi <<= 11;
				ct &= 0x7FF;
				ct |= gi;
				temp[count] = ct;
			}

			if ((temp[count] & 0x7) != NULL && (temp[count] & 0x7) != 0x1 && (temp[count] & 0x7) != 0x2)
				err->reportError(ERRTXT("Cell is no particle, nucleus or deformed structure "));

			c0 = ccurrent;
			count += 2;
			temp[count] = c0;
			temp[count + 1] = 1;
		}

	}

	if ((temp[count] & 0x7) != 0x2)
	{
		long ct = temp[count];
		long gi = ct >> 11;
		gi -= grains.defGrainsCount;

		if (gi < 0)
			err->reportError(ERRTXT("The grain index cannot be lower than 0"));

		gi <<= 11;
		ct &= 0x7FF;
		ct |= gi;
		temp[count] = ct;
	}

	long* finalArray = (long*)malloc((count + 2) * sizeof(long));
	if (!finalArray) err->reportError(ERRTXT("Cannot allocate more memory"));

	memcpy(finalArray, temp, sizeof(long) * (count + 2));
	free(temp);
	*size = (count + 2) / 2;
	return finalArray;
}

Real caHdl::nucleationProbability(Real rhoW) //###LB: Please check ths function  //as: not used
{
	//Real dsub = 0.0; //###LB: Carmen's subGrainDiameter2

	//if(gen.RVFlag == 1)
	//{
		//Do something. See Carmen's code
	//}
	//else
	//{
	//	dsub = parm.Cfac/sqrt(rho2W);
	//}

	//Real _nucProb =   dsub * val.subGrainDiameter * parm.Gam / ( parm.gam ) + ( 0.25 * val.subGrainDiameter* val.subGrainDiameter ); //###LB: Check parameters I could not find Gam or gam

	//Real nucProb = 1/_nucProb * parm.GBEfficiency; //###LB: This is what this function is supposed to deliver but please check the math.

	//parm.GBNucRest += (( N * grains.xSizeWanted * grains.ySizeWanted) - (long)(( N * grains.xSizeWanted * grains.ySizeWanted ))); //###LB: What is this? I think it is never used.

	//long N = (long) nucProb * grains.xSizeWanted * grains.ySizeWanted ; //###LB: Since grains are equally sized.

	Real nucProb = 0.01; //###LB: default for tests. Please change it once you have checked the equations.

	return nucProb;
}

cellP caHdl::rhoFacCP(float x)	//###LB: Changed	//rg: Not needed.	
{
	uint32_t ix = x * _RHOFACUNIT;
	QUICKASSERT(ix >= 0);
	QUICKASSERT(ix <= 255);
	//return (cellP) ( ix * 0x0100 + 1 ); 
	return (cellP)(ix * 0x00000008 + 1);
}

cellP caHdl::updateFacCP(float x, cellP c) // rg: Not needed.
{
	long ix = x * _RHOFACUNIT;
	QUICKASSERT(ix >= 0);
	QUICKASSERT(ix <= 255);
	long cp = ((long)c) & 0xFFFFF807;
	long newFac = ix * 0x0008;
	return (cellP)(cp | newFac);
}

grainP* caHdl::getGrainPos(long ix, long iy, long iz) // used during the original version when grains were rectangle. Not needed fo AM //as:not used
{
	QUICKASSERT(myNode() == parallelHdl->cellNode(ix, iy, iz)); //If exit; periodic

	ix = ix / grains.data.xCellCount;	// alles integer-division!    kein +0.5 !
	iy = iy / grains.data.yCellCount;
	iz = iz / grains.data.zCellCount;

	QUICKASSERT(ix == 0 && iy == 0 && iz == 0);

	return &grains.all[iz * grains.data.zFac + iy * grains.data.xPer + ix];
}

grainP caHdl::getOrigGrain(long ix, long iy, long iz)                                //Returns the corresponding grain for the selected coordinates
{
	QUICKASSERT(myNode() == parallelHdl->cellNode(ix, iy, iz)); //If exit; periodic

	ix = ix / grains.data.xCellCount;	// all integer division!		no +0.5!
	iy = iy / grains.data.yCellCount;
	iz = iz / grains.data.zCellCount;

	return grains.all[iz * grains.data.zFac + iy * grains.data.xPer + ix];
}

char caHdl::isParticle(cellP c)  //#as: checking for in sufficient melting before it was just return 0
{
	
	return 0 ; 
							

}

////////////////////
////////////////// porosity function

char caHdl::InsufficientMelting(cellP c) // declared in line 149 in caHdl.h//  inline long PorosityExists(cellP c) { return !((long)c & 0x00000001); } added in 51 grain.h
{											// it has to be void. i am checking 
	/*ofstream nomelt;
	nomelt.open("nomelt.txt");
	nomelt << "x\t\t\t y\t\t\t z\t\t\t T" << endl ;*/
	/*
	//cellP pd= (cell*)calloc((this->x0 * this->y0 * this->z0), sizeof(Real));
	//cellP porosity= NULL;/* = (cell*)calloc((this->x0 * this->y0 * this->z0), sizeof(Real));;
	vector<long> not_melted;
	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;
	long len = this->Tcells->T_timer.size;
	Real count = 0;
	Real number = 0;
	Real tpc = 0;

	for (uint32_t i = 0; i < (this->x0 * this->y0 * this->z0); i++)
	{
		pd[i].ix = 0;
		pd[i].iy = 0;
		pd[i].iz = 0;

		porosity[i].ix = 0;
		porosity[i].iy = 0;
		porosity[i].iz = 0;

	}*/

	
	//cout << "getNucTimeInactiveCell:" << timeAtSolidification << endl;
	//cout << "getCellNucleationTime:" << Tcells->getCellNucleationTime(c->dataHolder) << endl;
	//// for writing
	//porosity position;   //class porosity defined in grainca.h line 417
	long len = Tcells->T_timer.size - 1;
	Real maxTime = Tcells->T_timer.timeSeries[len];
	maxTime = maxTime * 2;
	Real tempt = Tcells->getCellNucleationTime(c->dataHolder);	//getNucTimeInactiveCell(c); //Tcells->getCellNucleationTime(c->dataHolder);
	if (tempt >= maxTime)
	{
		//cout << "insuff  " << c->ix << "\t" << c->iy << "\t" << c->iz << "\t" << endl;
		/*position.x = c->ix;
		position.y = c->iy;
		position.z = c->iz;
		position.total = nomelt.size();
		nomelt.push_back(position);
		position.total = nomelt.size();*/
		return 1;

	}
	else return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Output///////////////////////////////////


void caHdl::Output(void)
{
	if (rxOutput.now())
		RxFracOutput();

	OriOutput();
	MicrostructureOutput();
}

void caHdl::RxFracOutput(void)
{
	parm.rhomax = -1;	// unfug, müsste berechnet werden

	printf(" RXLocalFracTotal: %f time: %e\n", cells.data.totalRXfrac, this->time);
	printf("calc: %ld, %1.1lf MB\n", cells.data.count, cells.data.count * (sizeof(cell) / 1024.0 / 1024.4));
	printf("active:%ld\n", cells.data.count);
}

void caHdl::OriOutput(void)
{
	int mynode;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	int ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &ranks);

	Real rxFracArithMean;
	MPI_Allreduce(&(cells.data.totalRXfrac), &rxFracArithMean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	rxFracArithMean = rxFracArithMean / ranks;

	//###LB: See below the right order for this. If something may be empty or NULL you have to rule out first this possibility.
	//Otherwise it might lead to an addressing error. For instance, the code below.
	//You didn't see it because Linux can handle better some of such errors which is actually
	//wrong because it allows sloppy coding.

	if (!userDefinedOriFrequencies.empty() && rxFracArithMean >= userDefinedOriFrequencies.front()) {
		userDefinedOriFrequencies.pop_front();

		WriteOri();
		oriOutCount++;
	}
}

void caHdl::WriteOri(void)
{
	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	ostringstream fileNameStream;

	fileNameStream << "oribucket." << oriOutCount << ".ori";
	int fileNameLength = fileNameStream.str().size();
	char* fileName = new char[fileNameLength + 1];
	strcpy(fileName, fileNameStream.str().c_str());

	MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	int offset = 0;

	// write header
	if (myNode() == MASTER) {
		ostringstream headerString;
		headerString << "current simulation time t:/s:\t" << time << endl;
		headerString << "#oriIndex phi1 PHI phi2" << endl;
		headerString << "#Rank totalCellsInRank" << endl;
		headerString << "#oriIndex count" << endl;
		headerString << "TotalOris:\t" << grains.noris << endl;

		int headerStringLength = headerString.str().size();
		char* header = new char[headerStringLength + 1];
		strcpy(header, headerString.str().c_str());

		MPI_File_write_at(msFileHdl, offset, header, headerStringLength, MPI_CHAR, &msFileStatus);
		offset += headerStringLength;

		ostringstream totalOriString;
		for (int i = 0; i < grains.noris; i++) {
			totalOriString << i << "\t" << grains.orientations[i]->phi1 << "\t" << grains.orientations[i]->PHI << "\t" << grains.orientations[i]->phi2 << endl;
		}
		int totalOriStringLength = totalOriString.str().size();
		char* totalOri = new char[totalOriStringLength + 1];
		strcpy(totalOri, totalOriString.str().c_str());

		MPI_File_write_at(msFileHdl, offset, totalOri, totalOriStringLength, MPI_CHAR, &msFileStatus);
		offset += totalOriStringLength;
	}
	MPI_Bcast(&offset, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


	// write ori counts
	ostringstream nodeString;
	nodeString << "node:\t" << myNode() << "\t" << xPer * yPer * zPer << endl;
	for (int i = 0; i < grains.noris; i++) {
		nodeString << i << "\t" << grains.orientations[i]->count;
		if (i != grains.noris - 1) nodeString << endl;
	}

	// looks weird but has to be done to ensure uniform string length for MPI IO
	int nodeStringLength = nodeString.str().size();
	if (nodeStringLength < (28 + 21 * grains.noris - 1)) { // size of the block, fill with spaces
		int spaceCount = 28 + 22 * grains.noris - nodeStringLength - 1;
		for (int i = 0; i < spaceCount; i++) {
			nodeString << " ";
		}
	}
	nodeString << endl;

	char* nodeStr = new char[(28 + 22 * grains.noris) + 1];
	strcpy(nodeStr, nodeString.str().c_str());

	offset += myNode() * (28 + 22 * grains.noris);
	MPI_File_write_at(msFileHdl, offset, nodeStr, 28 + 22 * grains.noris, MPI_CHAR, &msFileStatus);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&msFileHdl);
	MPI_Barrier(MPI_COMM_WORLD);
}

void caHdl::MicrostructureOutput(void)  // use these functions as per requirements

{
	int mynode;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	int ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &ranks);

	// User defined Frequencies
	if (userDefinedMSOutput) {
		Real rxFracArithMean;
		MPI_Allreduce(&(cells.data.totalRXfrac), &rxFracArithMean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//collect localized recrystallized fraction WITH CONSIDERATION OF PARTIALLY INFECTED CELLS AND RECRYSTALLIZED FRACTION WHICH
		//WHICH HAS BEEN GIVEN RECRYSTALLIZED CELLS IN ADVANCE BECAUSE OF THE JUMP CONSTRUCTION IN THE MAIN LOOP into rxFracArithMean
		//combine from all processes and distributed back to all processes
		//SHOULD BE PARTICULARLY ADDRESS WHEN PARTICLES ARE CONSIDERED BECAUSE ONLY THE NON RECRYSTALLIZABLE VOLUME IS ACCOUNTED FOR####
		rxFracArithMean = rxFracArithMean / ranks;

		// total recrystallisation is bigger or equal to the smallest user defined amount, if list is not empty and if plot shall be done!
	   // if ( !userDefinedFrequencies.empty() && rxFracArithMean >= userDefinedFrequencies.front() && generateMSBinary ) {
		if (!userDefinedFrequencies.empty() && time >= userDefinedFrequencies.front() && generateMSBinary) {// was written for sepi, output at particular time intervals
			userDefinedFrequencies.pop_front();

			/*if (this->MSBinaryPlotMode == PLOTINFPROGRESS)
				WriteMicrostructureInfGreyscale();	//GREYSCALE MAP
			else
				WriteMicrostructurePRNGRGB();*/													//ORIINDEX MAP

				//write_xsection_coloring_ipfz();
			//	write_voxeldata_coloring_ipfz();
				//write_voxeldata_oris();
			//draw_pole_figure();
			//writingVoxelOrientationID_SEPI();
			//write_state_SEPI();

			outCount++;
		}
	}
	// Fixed Frequencies
	else {
		// not yet implemented
	}
}

//###Check correctness for the set file view snippet
//File_seek operation can shift not more than 2^32 Bytes per operation appr. 4.2billion bytes
//858Mio IOCellObjects
// Hoefter 08102012 (takes in account the possibility of differently sized nodes)
void caHdl::WriteMicrostructureInfGreyscale(void) //to allow plotting the infection state of each recrystallized cell ranging from black (0) to white (255)
{
	int mynode;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);       //who am I?
	int ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &ranks);         //how many are out there?

	// Write Microstructure using Parallel IO
	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O
	ostringstream fileNameStream;
	fileNameStream << "microstructure_transfprogress." << outCount;
	int fileNameLength = fileNameStream.str().size();
	char* fileName = new char[fileNameLength + 1];
	strcpy(fileName, fileNameStream.str().c_str());

	// open the file in create and write-only mode
	MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);
	int totalOffset = 0;

	// ori indices are not necessary
	if (mynode == MASTER) {
		int nodeCount = ranks;
		//write number of nodes
		MPI_File_write_at(msFileHdl, totalOffset, &nodeCount, 1, MPI_INT, &msFileStatus);
		totalOffset += 4; // + int
		//write cellsize micron/1cell
		double cellSize = cells.data.size;
		MPI_File_write_at(msFileHdl, totalOffset, &cellSize, 1, MPI_DOUBLE, &msFileStatus);
		totalOffset += 8; // + double
	}

	// Broadcast totalOffset the position after which the node-independent header information have ended in memory
	MPI_Bcast(&totalOffset, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	totalOffset += mynode * 28; // 28 = sizeof(MPI_IO_Node_Type), what if there are so many nodes?###

	// fill MPI_IO_Node object with my node-specific data
	MPI_IO_Node nodeLine;
	nodeLine.MPIRank = mynode;
	nodeLine.x0 = this->x0;
	nodeLine.y0 = this->y0;
	nodeLine.z0 = this->z0;
	nodeLine.xmax = this->x0 + this->xPer;
	nodeLine.ymax = this->y0 + this->yPer;
	nodeLine.zmax = this->z0 + this->zPer;

	// write to file
	MPI_File_write_at(msFileHdl, totalOffset, &nodeLine, 1, MPI_IO_Node_Type, &msFileStatus);
	// "sync" offset
	totalOffset += (ranks - mynode) * 28;

	// create xyzOffset array (one for every rank)
	int* xyzOffset = new int[ranks];
	xyzOffset[mynode] = this->xPer * this->yPer * this->zPer; //###non portable only 4.28 Billions cells per Node

	// every node broadcasts its value, performance bottleneck...
	for (int i = 0; i < ranks; i++) {
		MPI_Bcast(xyzOffset + i, 1, MPI_INT, i, MPI_COMM_WORLD);
	}

	// seek pointer to start writing the my binary block with my cells out
	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);
	for (int i = 0; i < mynode; i++) {
		// for every byte of an MPI_IO_Cell, because this can easily exceed the limits of an int
		for (int j = 0; j < 5; j++) {
			MPI_File_seek(msFileHdl, xyzOffset[i], MPI_SEEK_CUR); // instead of this->xPer*this->yPer*this->zPer, access the xyzOffset of the i'th rank
		}
	}

	// iterate through all xy-Planes
	for (int z = 0; z < this->zPer; z++) {
		MPI_IO_CellInfState* cellPlane = new MPI_IO_CellInfState[this->yPer * this->xPer]; //this could also be saved...

		for (int y = 0; y < this->yPer; y++) {
			for (int x = 0; x < this->xPer; x++) {
				int cpIndex = y * this->xPer + x;
				cellP currentCell = getCellPer(this->x0 + x, this->y0 + y, this->z0 + z);

				//###categorize a cell, everything which is not zero is true, maybe here it is possible to save some if clauses if the order is changed...
				if (currentCell != NULL) { //cell has been created but assigned either a particle or a recrystallizing grain
					cellPlane[cpIndex].InfProgress = (float)(255 * (1.0 - currentCell->rxFrac)); //##caveat: also cell which have been infected <= 3.92E-3 are undistinguishable
					//###there should be a routine which checks in how far the unsigned interval is not exceeded!
					//THE FATAL EXCEPTION cellPlane[cpIndex].InfProgress < 0 is checkec by mstoPNG
					cellPlane[cpIndex].state = 1; //cell is representing an infected cell
				}
				else { //cell has not been created as it is still deformed
					cellPlane[cpIndex].InfProgress = 0.0;
					cellPlane[cpIndex].state = 0; //cell is representing a part of a uninfected deformed grain
				}
			}
		}

		// write the memory block, increases the file-pointer implicitly, blocking, non-collectively 
		MPI_File_write(msFileHdl, cellPlane, this->xPer * this->yPer, MPI_IO_CellInfState_Type, &msFileStatus);

		delete[] cellPlane; // cleanup
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&msFileHdl);

}

//###Check correctness for the set file view snippet
//File_seek operation can shift not more than 2^32 Bytes per operation appr. 4.2billion bytes
//858Mio IOCellObjects
// Hoefter 08102012 (takes in account the possibility of differently sized nodes)
void caHdl::WriteMicrostructurePRNGRGB(void) //to allow plotting the orientation indices of each cell
{
	int mynode = rank;
	//MPI_Comm_rank( MPI_COMM_WORLD, &mynode );       //who am I?
	int ranks = totalNodes;
	//MPI_Comm_size( MPI_COMM_WORLD, &ranks);         //how many are out there?

	// Write Microstructure using Parallel IO
	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O
	ostringstream fileNameStream;
	fileNameStream << "microstructure_oriindices." << outCount;
	int fileNameLength = fileNameStream.str().size();
	char* fileName = new char[fileNameLength + 1];
	strcpy(fileName, fileNameStream.str().c_str());

	// open the file in create and write-only mode
	MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);
	int totalOffset = 0;

	// only Master needs to write oriindices
	if (mynode == MASTER) {
		int oriCount = grains.noris;
		int nodeCount = ranks;

		MPI_IO_OriIndex* indices = new MPI_IO_OriIndex[oriCount]; // allocate OriIndex array
		for (int i = 0; i < oriCount; i++) {
			indices[i].OriIndex = i;
			indices[i].phi1 = (float)grains.orientations[i]->phi1;
			indices[i].PHI = (float)grains.orientations[i]->PHI;
			indices[i].phi2 = (float)grains.orientations[i]->phi2;
		}

		//header of the microstructure binary: int<nodecount>int<oricount> (double<cellSize>)
		MPI_File_write_at(msFileHdl, totalOffset, &nodeCount, 1, MPI_INT, &msFileStatus);
		totalOffset += 4; // + int
		MPI_File_write_at(msFileHdl, totalOffset, &oriCount, 1, MPI_INT, &msFileStatus);
		totalOffset += 4; // + int
		//if (INCLUDECELLSIZE) {  // Required for mstoGrainSize (and new versions of mstoPNG) provides dimension of the cubic cell to allow the digital line interception method
		double cellSize = cells.data.size;
		MPI_File_write_at(msFileHdl, totalOffset, &cellSize, 1, MPI_DOUBLE, &msFileStatus);
		totalOffset += 8; // + double
		//}

		//write memory block with contains mpiiooriindices ...MPIIOOriIndex<>*oriCount
		MPI_File_write_at(msFileHdl, totalOffset, indices, oriCount, MPI_IO_OriIndex_Type, &msFileStatus);
		totalOffset += oriCount * 16; // 16 = sizeof(MPI_IO_OriIndex)

		delete[] indices; // cleanup
	}

	// Broadcast totalOffset the position after which the node-independent header information have ended in memory
	MPI_Bcast(&totalOffset, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	totalOffset += mynode * 28; // 28 = sizeof(MPI_IO_Node_Type), what if there are so many nodes?###

	// fill MPI_IO_Node object with my node-specific data
	MPI_IO_Node nodeLine;
	nodeLine.MPIRank = mynode;
	nodeLine.x0 = this->x0;
	nodeLine.y0 = this->y0;
	nodeLine.z0 = this->z0;
	nodeLine.xmax = this->x0 + this->xPer;
	nodeLine.ymax = this->y0 + this->yPer;
	nodeLine.zmax = this->z0 + this->zPer;

	// write to file
	MPI_File_write_at(msFileHdl, totalOffset, &nodeLine, 1, MPI_IO_Node_Type, &msFileStatus);
	// "sync" offset
	totalOffset += (ranks - mynode) * 28;

	// create xyzOffset array (one for every rank)
	int* xyzOffset = new int[ranks];
	xyzOffset[mynode] = this->xPer * this->yPer * this->zPer; //###non portable only 4.28 Billions cells per Node

	// every node broadcasts its value, performance bottleneck...
	for (int i = 0; i < ranks; i++) {
		MPI_Bcast(xyzOffset + i, 1, MPI_INT, i, MPI_COMM_WORLD);
	}

	// seek pointer to start writing the my binary block with my cells out
	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);
	for (int i = 0; i < mynode; i++) {
		// for every byte of an MPI_IO_Cell, because this can easily exceed the limits of an int
		for (int j = 0; j < 5; j++) {
			MPI_File_seek(msFileHdl, xyzOffset[i], MPI_SEEK_CUR); // instead of this->xPer*this->yPer*this->zPer, access the xyzOffset of the i'th rank
		}
	}

	// iterate through all xy-Planes
	for (int z = 0; z < this->zPer; z++) {
		MPI_IO_Cell* cellPlane = new MPI_IO_Cell[this->yPer * this->xPer]; //this could also be saved...

		for (int y = 0; y < this->yPer; y++) {
			for (int x = 0; x < this->xPer; x++) {
				int cpIndex = y * this->xPer + x;
				cellP currentCell = getCellPer(this->x0 + x, this->y0 + y, this->z0 + z);

				//###categorize a cell, everything which is not zero is true, maybe here it is possible to save some if clauses if the order is changed...
				if (currentCell != NULL) { //cell has been created but assigned either a particle or a recrystallizing grain
					grainP currentGrain = grainOfCell(currentCell);
					cellPlane[cpIndex].OriIndex = currentGrain->oriIndex;
					cellPlane[cpIndex].state = 1; //cell is representing a infectedOr recrystallized grain###
				}
				else {
					cellPlane[cpIndex].OriIndex = getOrigGrain(this->x0 + x, this->y0 + y, this->z0 + z)->oriIndex;
					cellPlane[cpIndex].state = 0; //cell is representing a part of a uninfected deformed grain
				}
			}
		}

		// write the memory block, increases the file-pointer implicitly, blocking, non-collectively 
		MPI_File_write(msFileHdl, cellPlane, this->xPer * this->yPer, MPI_IO_Cell_Type, &msFileStatus);

		delete[] cellPlane; // cleanup
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&msFileHdl);

}

void caHdl::write_xsection_coloring_ipfz(void)
{
	int myRank = myNode();
	long counter = 0;
	unsigned char* img_rgba = NULL;
	ostringstream fname;

	//gen.backgroundMode = BLACK;  //It must be read from input file.

	uint32_t imgx = this->box.xPer;    //myCAGeometry.nboxedge_rd;
	uint32_t imgy = this->box.yPer;    //myCAGeometry.nboxedge_nd;
	uint32_t imgxy = imgx * imgy;

	uint32_t imgz = (uint32_t)(DEFAULT_ZSECTIONING_ZPOS * ((double)this->box.zPer));

	if (myRank == MASTER)
	{
		img_rgba = (unsigned char*)calloc(imgx * imgy * 4, sizeof(unsigned char));

		fname << "SCORE_X.";
		uint32_t XX = (cells.data.totalRXfrac * 1000.0);
		if (XX < 10)					fname << "000" << XX;
		if (XX >= 10 && XX < 100)		fname << "00" << XX;
		if (XX >= 100 && XX < 1000)	fname << "0" << XX;
		if (XX >= 1000)				fname << XX;
		fname << ".IPFZ.png";

		cout << "Rendering a section__" << fname.str().c_str() << "__ at " << imgx << ";" << imgy << ";" << imgz << endl;

	}

	long x = (long)(this->box.zPer * 0.5);

	for (long y = 0; y < this->box.yPer; y++) {
		for (long z = 0; z < this->box.zPer; z++) {

			MPI_Barrier(MPI_COMM_WORLD);

			int cpIndex = y * this->xPer + x;

			cellP currentCell = NULL;

			unsigned char colors[4] = { WHITE_RGB,WHITE_RGB,WHITE_RGB,255 }; //default is white

			int cellsNode = parallelHdl->cellNode(x, y, z);

			if (myRank == cellsNode)
			{
				cellP* cp = getCellPerP(x, y, z);
				cellP c = *cp;

				if (!cellExists(c)) //Cell belongs to deformed grain
				{
					if (gen.backgroundMode == DEFMODUS)
					{
						grainP g = getOrigGrain(c);
						long localOriIndex = g->oriIndex;
						oriRepresentationP localOri = grains.orientations[localOriIndex];
						colors[RED] = localOri->RGBA[RED];
						colors[GREEN] = localOri->RGBA[GREEN];
						colors[BLUE] = localOri->RGBA[BLUE];
						colors[ALPHA] = localOri->RGBA[ALPHA];
					}
					else if (gen.backgroundMode == BLACK)
						memset(colors, 0, 3 * sizeof(unsigned char));
				}
				else if (isParticle(c))
				{
					colors[RED] = BLACK;
					colors[GREEN] = BLACK;
					colors[BLUE] = BLACK;
					colors[ALPHA] = BLACK;
				}
				else
				{
					grainP g = grainOfCell(c);
					long localOriIndex = g->oriIndex;
					oriRepresentationP localOri = grains.orientations[localOriIndex];

					if (isCellRecrystallized(c))
					{
						colors[RED] = localOri->RGBA[RED];
						colors[GREEN] = localOri->RGBA[GREEN];
						colors[BLUE] = localOri->RGBA[BLUE];
						colors[ALPHA] = localOri->RGBA[ALPHA];
					}
					else
					{
						colors[RED] = (unsigned char)(0.5 * ((float)localOri->RGBA[RED] + (float)WHITE_RGB));
						colors[GREEN] = (unsigned char)(0.5 * ((float)localOri->RGBA[GREEN] + (float)WHITE_RGB));
						colors[BLUE] = (unsigned char)(0.5 * ((float)localOri->RGBA[BLUE] + (float)WHITE_RGB));
						colors[ALPHA] = localOri->RGBA[ALPHA];
					}
				}

				if (myRank != MASTER)
				{
					MPI_Send(&colors, 4, MPI_UNSIGNED_CHAR, MASTER, myNode(), MPI_COMM_WORLD);
				}
				else if (myRank == MASTER)
				{
					img_rgba[counter * 4 + RED] = colors[RED];
					img_rgba[counter * 4 + GREEN] = colors[GREEN];
					img_rgba[counter * 4 + BLUE] = colors[BLUE];
					img_rgba[counter * 4 + ALPHA] = colors[ALPHA];
					counter++;
				}

			}
			else if (myRank == MASTER)
			{
				unsigned char rcvData[4] = { 0 };
				MPI_Status status;
				MPI_Recv(&rcvData, 4, MPI_UNSIGNED_CHAR, cellsNode, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				img_rgba[counter * 4 + RED] = rcvData[RED];
				img_rgba[counter * 4 + GREEN] = rcvData[GREEN];
				img_rgba[counter * 4 + BLUE] = rcvData[BLUE];
				img_rgba[counter * 4 + ALPHA] = rcvData[ALPHA];
				counter++;
			}
		}
	}

	if (myRank == MASTER)	lodepng::encode(fname.str().c_str(), img_rgba, imgx, imgy);

	free(img_rgba);
}

void caHdl::write_voxeldata_coloring_ipfz(void)
{
	//TARGET FORMAT IS INT FOR PARAVIEW OR AVIZO ONLY RECRYSTALLIZED GRAINS (rxFrac < 1.0) ARE colored white
	//CELL_IS_PARTICLE	0
	//CELL_IS_INFECTED	1
	//DEFORMED 			1+mydefgid
	//RECRYSTALLIZED	1+mydefgid.size()+myrxgid
	//<MPI_INT>uint

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O

	stringstream msFileName;
	msFileName << "SCORE."  "time." << this->time;
	msFileName << ".MS3DIPFZ.raw";

	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength + 1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl); // changing MPI_COMM_SELF to  MPI_COMM_WORLD

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized

	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;

	unsigned char* xRowData = (unsigned char*)calloc(RGB * (nx - x0), sizeof(unsigned char));

	long long totalOffset = 0;

	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);

	MPI_Barrier(MPI_COMM_WORLD);

	for (uint32_t z = z0; z < nz; z++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * RGB;
		for (uint32_t y = y0; y < ny; y++)
		{
			long counterx = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			for (uint32_t x = x0; x < nx; x++)
			{
				cellP currentCell = NULL;
				unsigned char rgb[3] = { WHITE_RGB, WHITE_RGB, WHITE_RGB };

				cellP* cp = getCellPerP(x, y, z);
				cellP c = *cp;
				oriRepresentation oriDummy;

				if (cellExists(c))
				{
					if (isParticle(c) || InsufficientMelting(c))
					{
						rgb[RED] = BLACK;
						rgb[GREEN] = BLACK;
						rgb[BLUE] = BLACK;
					}
					else
					{
						grainP g = grainOfCell(c);
						long localOriIndex = g->oriIndex;
						oriRepresentationP localOri = grains.orientations[localOriIndex];
						oriDummy.copyOriFrom(localOri);
						rgb[RED] = (unsigned char)(0.5 * ((float)WHITE_RGB + (float)oriDummy.RGBA[RED]));
						rgb[GREEN] = (unsigned char)(0.5 * ((float)WHITE_RGB + (float)oriDummy.RGBA[GREEN]));
						rgb[BLUE] = (unsigned char)(0.5 * ((float)WHITE_RGB + (float)oriDummy.RGBA[BLUE]));

						if (isCellRecrystallized(c))
						{
							rgb[RED] = oriDummy.RGBA[RED];
							rgb[GREEN] = oriDummy.RGBA[GREEN];
							rgb[BLUE] = oriDummy.RGBA[BLUE];
						}
					}
				}

				xRowData[counterx * RGB + RED] = rgb[RED];
				xRowData[counterx * RGB + GREEN] = rgb[GREEN];
				xRowData[counterx * RGB + BLUE] = rgb[BLUE];
				counterx++;
			}
			MPI_File_write_at(msFileHdl, totalOffset, xRowData, counterx * RGB, MPI_UNSIGNED_CHAR, &msFileStatus);
			totalOffset += Nx * RGB;
		}
	}

	delete[] CmsFileName;
	free(xRowData);

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}

void caHdl::write_voxeldata_oris(void)
{
	//TARGET FORMAT IS INT FOR PARAVIEW OR AVIZO ONLY RECRYSTALLIZED GRAINS (rxFrac < 1.0) ARE colored white
//CELL_IS_PARTICLE	0
//CELL_IS_INFECTED	1
//DEFORMED 			1+mydefgid
//RECRYSTALLIZED	1+mydefgid.size()+myrxgid
//<MPI_INT>uint

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O

	stringstream msFileName;
	msFileName << "SCORE."  "time." << this->time;
	msFileName << ".EulerOris.raw";

	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength + 1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized

	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;

	unsigned char* xRowData = NULL;
	xRowData = new unsigned char[(nx - x0) * 12]; // 3 * (4 bytes) numbers 

	long long totalOffset = 0;

	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);

	MPI_Barrier(MPI_COMM_WORLD);

	unsigned char bytes[12] = { 0 };  // 3 numbers each with 4 bytes = 12 bytes

	for (uint32_t z = z0; z < nz; z++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * 12;
		for (uint32_t y = y0; y < ny; y++)
		{
			long counterx = 0;
			MPI_Barrier(MPI_COMM_WORLD);
			for (uint32_t x = x0; x < nx; x++)
			{
				cellP currentCell = NULL;
				uint32_t rgb[3] = { (uint32_t)0, (uint32_t)0, (uint32_t)0 };

				cellP* cp = getCellPerP(x, y, z);
				cellP c = *cp;
				oriRepresentation oriDummy;

				if (cellExists(c))
				{
					if (InsufficientMelting(c)|| isParticle(c))
					{
						rgb[RED] = (uint32_t)BLACK;
						rgb[GREEN] = (uint32_t)BLACK;
						rgb[BLUE] = (uint32_t)BLACK;
					}
					else
					{
						grainP g = grainOfCell(c);
						long localOriIndex = g->oriIndex;
						oriRepresentationP localOri = grains.orientations[localOriIndex];
						oriDummy.copyOriFrom(localOri);

						rgb[RED] = (uint32_t)(returnPositiveAngle(oriDummy.phi1) * 10000);
						rgb[GREEN] = (uint32_t)(returnPositiveAngle(oriDummy.PHI) * 10000);
						rgb[BLUE] = (uint32_t)(returnPositiveAngle(oriDummy.phi2) * 10000);

						if (isCellRecrystallized(c))
						{
							rgb[RED] = (uint32_t)(returnPositiveAngle(oriDummy.phi1) * 10000);
							rgb[GREEN] = (uint32_t)(returnPositiveAngle(oriDummy.PHI) * 10000);
							rgb[BLUE] = (uint32_t)(returnPositiveAngle(oriDummy.phi2) * 10000);
						}
					}
				}

				for (int32_t i = 0; i < 3; i++)
					for (int32_t j = 3; j >= 0; j--)
					{
						//bytes[(3 - j) + 4 * i] = (rgb[i] >> (8 * j)) & 0xFF;// storing data as big endian form.
						bytes[j + 4 * i] = (rgb[i] >> (8 * j)) & 0xFF;// storing data as little endian form.
					}

				memcpy(&xRowData[counterx * 12], bytes, 12 * sizeof(unsigned char));

				counterx++;
			}

			MPI_File_write_at(msFileHdl, totalOffset, xRowData, 12 * counterx, MPI_UNSIGNED_CHAR, &msFileStatus);
			totalOffset += ((long long)Nx * 12);
		}
	}

	delete[] CmsFileName;
	delete[] xRowData;

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// WRITE FUNCTIONS TO PRINT SEPI DATA.

void caHdl::orientations_SEPI(void) {
	// write orientation index and the theta values.
	ofstream orientaionsFile("Orientaions_SEPI", ios::out);
	orientaionsFile << "index" << "\t" << "phi_1" << "\t" << "PHI" << "\t" << "phi_2" << "\t" << endl;

	if (!orientaionsFile) {
		cout << "File could not be created" << endl;
	}
	for (int i = 0; i < this->grains.noris; i++) {

		double local_phi_1 = this->grains.orientations[i]->phi1;
		double local_phi = this->grains.orientations[i]->PHI;
		double local_phi_2 = this->grains.orientations[i]->phi2;

		orientaionsFile << i << "\t" << local_phi_1 << "\t" << local_phi << "\t" << local_phi_2 << endl;
	}

	orientaionsFile.close();
}

void caHdl::writingVoxelOrientationID_SEPI()
{

	//TARGET FORMAT IS INT FOR PARAVIEW OR AVIZO ONLY RECRYSTALLIZED GRAINS (rxFrac < 1.0) ARE colored white
	//CELL_IS_PARTICLE	0
	//CELL_IS_INFECTED	1
	//DEFORMED 			1+mydefgid
	//RECRYSTALLIZED	1+mydefgid.size()+myrxgid
	//<MPI_INT>uint

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O

	stringstream msFileName;
	msFileName << "SCORE."  "time." << this->time;
	msFileName << ".VoxelOrientaionIndex_SEPI.raw";

	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength + 1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized

	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;

	unsigned char* xRowData = NULL;
	xRowData = new unsigned char[(nx - x0) * 4]; // 4 bytes to represent orientation index 

	long long totalOffset = 0;

	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);

	MPI_Barrier(MPI_COMM_WORLD);

	unsigned char bytes[4] = { 0 };

	for (uint32_t z = z0; z < nz; z++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * 4;
		for (uint32_t y = y0; y < ny; y++)
		{
			long counterx = 0;
			MPI_Barrier(MPI_COMM_WORLD);

			for (uint32_t x = x0; x < nx; x++)
			{
				cellP currentCell = NULL;
				uint32_t orientationIndex = 0;

				cellP* cp = getCellPerP(x, y, z);
				cellP c = *cp;

				if (cellExists(c))
				{
					grainP g = grainOfCell(c);
					orientationIndex = g->oriIndex;
				}

				for (int32_t j = 3; j >= 0; j--)
				{
					//bytes[(3 - j) +] = (orientationIndex >> (8 * j)) & 0xFF;// storing data as big endian form.
					bytes[j] = (orientationIndex >> (8 * j)) & 0xFF;// storing data as little endian form.
				}

				memcpy(&xRowData[counterx * 4], bytes, 4 * sizeof(unsigned char));

				counterx++;
			}

			MPI_File_write_at(msFileHdl, totalOffset, xRowData, 4 * counterx, MPI_UNSIGNED_CHAR, &msFileStatus);
			totalOffset += ((long long)Nx * 4);
		}
	}

	delete[] CmsFileName;
	delete[] xRowData;

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF

}

void caHdl::write_state_SEPI(void) //  rg: this function was added for sepi data.
{
#define REGDATSIZ (5)
	// #define MELTINGPOINT (1375) ## as now reading from input file

#define CELL_IS_UNKNOWN 0
#define	CELL_IS_POWDER 1;
#define CELL_IS_LIQUID  2;
#define CELL_IS_PARTIALLY_CRYSTALLIZED 3;
#define CELL_IS_CRYSTALLIZED 4;
#define CELL_IS_PARTICLE 5;
	Real MELTINGPOINT = gen.Tmelt;

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O

	stringstream msFileName;
	msFileName << "SCORE."  "time." << this->time;
	//msFileName << ".MS3DIPFZ.raw";
	msFileName << ".MS3D_STATE_SEPI.raw";

	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength + 1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized

	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;

	//unsigned char* xRowData = (unsigned char*)calloc(RGB * (nx - x0), sizeof(unsigned char));
	unsigned char* xRowData = (unsigned char*)calloc((nx - x0), sizeof(unsigned char));

	long long totalOffset = 0;

	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);

	MPI_Barrier(MPI_COMM_WORLD);

	long CELL_IS_POWDER_count = 0;
	long CELL_IS_LIQUID_count = 0;
	long CELL_IS_PARTIALLY_CRYSTALLIZED_count = 0;
	long CELL_IS_CRYSTALLIZED_count = 0;
	long CELL_IS_PARTICLE_count = 0;
	long CELL_IS_UNKNOWN_count = 125000000;

	for (uint32_t z = z0; z < nz; z++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		//totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * RGB;
		totalOffset = (x0 + z * Nx * Ny + y0 * Nx);

		for (uint32_t y = y0; y < ny; y++)
		{
			long counterx = 0;
			MPI_Barrier(MPI_COMM_WORLD);

			for (uint32_t x = x0; x < nx; x++)
			{
				cellP currentCell = NULL;
				//unsigned char rgb[3] = { WHITE_RGB, WHITE_RGB, WHITE_RGB };

				//uint8_t state = CELL_IS_UNKNOWN;
				unsigned char state = CELL_IS_UNKNOWN;

				cellP* cp = getCellPerP(x, y, z);
				cellP c = *cp;
				//oriRepresentation oriDummy;

				if (!cellExists(c))
				{

					tempCellP tc = getTCellofCell(c);
					//tempCellP tc = Tcells->PNBucket.TCparallelNode[i];

					short x0 = tc->mySpace.xleft; // these are the beginning and the end of the temperature_cells in all direction.
					short xf = tc->mySpace.xright;
					short y0 = tc->mySpace.yleft;
					short yf = tc->mySpace.yright;
					short z0 = tc->mySpace.zleft;
					short zf = tc->mySpace.zright;

					//Real ix = x;
					//Real iy = y;
					//Real iz = z;

					Real ix = (x0 + xf) / 2;
					Real iy = (y0 + yf) / 2;
					Real iz = (z0 + zf) / 2;

					Real xd = (Real)(ix - x0) / (Real)(xf - x0);
					Real yd = (Real)(iy - y0) / (Real)(yf - y0);
					Real zd = (Real)(iz - z0) / (Real)(zf - z0);

					Real currentNodesTemperature[8] = { 0 };
					TNodeP TCellNodes[8] = { tc->data.vertexTemp[c000], tc->data.vertexTemp[c001], tc->data.vertexTemp[c010], tc->data.vertexTemp[c011],
					 tc->data.vertexTemp[c100], tc->data.vertexTemp[c101], tc->data.vertexTemp[c110], tc->data.vertexTemp[c111] };

					Real t[5] = { 0 };
					Real T[5] = { 0 };
					Real m = 0; //slope
					Real c0 = 0; //intercept
					Real R2 = 0; // R-squared
					Real meltTime = NO_MELT;

					for (long j = 0; j < REGDATSIZ; j++)
					{
						long timeKey = Tcells->T_timer.key + j;
						long nextKey = timeKey + 1;

						if (timeKey >= Tcells->T_timer.size - 1)
						{
							timeKey = Tcells->T_timer.size - 1;
							nextKey = timeKey;
						}

						double previous_time = Tcells->T_timer.timeSeries[timeKey];
						double current_time = this->time;
						double next_time = Tcells->T_timer.timeSeries[nextKey];

						t[j] = previous_time;

						for (int i = c000; i <= c111; i++)
							currentNodesTemperature[i] = TCellNodes[i]->interpolateTemperatureTime(previous_time, current_time, next_time, timeKey, nextKey);

						Real T00 = currentNodesTemperature[c000] * (1 - xd) + currentNodesTemperature[c100] * xd;
						Real T01 = currentNodesTemperature[c001] * (1 - xd) + currentNodesTemperature[c101] * xd;
						Real T10 = currentNodesTemperature[c010] * (1 - xd) + currentNodesTemperature[c110] * xd;
						Real T11 = currentNodesTemperature[c011] * (1 - xd) + currentNodesTemperature[c111] * xd;
						Real T0 = T00 * (1 - yd) + T10 * yd;
						Real T1 = T01 * (1 - yd) + T11 * yd;
						T[j] = T0 * (1 - zd) + T1 * zd;

					}

					if (T[0] >= MELTINGPOINT && T[REGDATSIZ] >= MELTINGPOINT)
					{
						state = CELL_IS_LIQUID;
						CELL_IS_LIQUID_count++;
						CELL_IS_UNKNOWN_count--;
					}

					else {
						linearRegression(t, T, REGDATSIZ, &m, &c0, &R2);
						Real Tcrit = m * t[REGDATSIZ - 1] + c0;
						if (Tcrit < MELTINGPOINT)
						{
							state = CELL_IS_POWDER;
							CELL_IS_POWDER_count++;
							CELL_IS_UNKNOWN_count--;
						}

						else
						{
							state = CELL_IS_LIQUID;
							CELL_IS_LIQUID_count++;
							CELL_IS_UNKNOWN_count--;

						}

					}
				}

				if (cellExists(c))
				{
					//cout << "exist" << endl;
					if (isCellRecrystallized(c))
					{
						state = CELL_IS_CRYSTALLIZED;
						CELL_IS_CRYSTALLIZED_count++;
						CELL_IS_UNKNOWN_count--;
					}

					else if (c->rxFrac == 0) {
						state = CELL_IS_POWDER;
						CELL_IS_POWDER_count++;
						CELL_IS_UNKNOWN_count--;
					}

					else
					{
						state = CELL_IS_PARTIALLY_CRYSTALLIZED;
						CELL_IS_PARTIALLY_CRYSTALLIZED_count++;
						CELL_IS_UNKNOWN_count--;
					}
				}

				xRowData[counterx] = state;
				counterx++;
			}
			//MPI_File_write_at(msFileHdl, totalOffset, xRowData, counterx * RGB, MPI_UNSIGNED_CHAR, &msFileStatus);
			//totalOffset += Nx * RGB;
			MPI_File_write_at(msFileHdl, totalOffset, xRowData, counterx, MPI_UNSIGNED_CHAR, &msFileStatus);
			totalOffset += Nx;
		}
	}

	delete[] CmsFileName;
	free(xRowData);

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF

	cout << "KNOWN" << "\t" << CELL_IS_UNKNOWN_count << endl;
	cout << "POWDER" << "\t" << CELL_IS_POWDER_count << endl;
	cout << "LIQUID" << "\t" << CELL_IS_LIQUID_count << endl;
	cout << "PARTIALLY_CRYSTALLIZED" << "\t" << CELL_IS_PARTIALLY_CRYSTALLIZED_count << endl;
	cout << "CRYSTALLIZED" << "\t" << CELL_IS_CRYSTALLIZED_count << endl;
	cout << "PARTICLE" << "\t" << CELL_IS_PARTICLE_count << endl;

	cout << " total_count" << CELL_IS_UNKNOWN_count + CELL_IS_POWDER_count + CELL_IS_LIQUID_count + CELL_IS_PARTIALLY_CRYSTALLIZED_count + CELL_IS_CRYSTALLIZED_count + CELL_IS_PARTICLE_count << endl;
}


void caHdl::write_final_state_SEPI(void) //  rg: this function was added for sepi data.
{
#define REGDATSIZ (5)
	// #define MELTINGPOINT (1350)  #as : Now reading from input file

	float MELTINGPOINT = gen.Tmelt; // #as : Now reading from input file


#define CELL_IS_UNKNOWN 0
#define	CELL_IS_POWDER 1;
#define CELL_IS_LIQUID  2;
#define CELL_IS_PARTIALLY_CRYSTALLIZED 3;
#define CELL_IS_CRYSTALLIZED 4;
#define CELL_IS_PARTICLE 5;
#define CELL_IS_UnmeltedPowder 6;

	MPI_File msFileHdl;
	MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O

	stringstream msFileName;
	msFileName << "SCORE."  "time." << this->time;
	//msFileName << ".MS3DIPFZ.raw";
	msFileName << ".MS3D_STATE_SEPI.raw";

	int msFileNameLength = msFileName.str().size();
	char* CmsFileName = new char[msFileNameLength + 1];
	strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized

	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;

	//unsigned char* xRowData = (unsigned char*)calloc(RGB * (nx - x0), sizeof(unsigned char));
	unsigned char* xRowData = (unsigned char*)calloc((nx - x0), sizeof(unsigned char));

	long long totalOffset = 0;

	MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);

	MPI_Barrier(MPI_COMM_WORLD);
	long CELL_IS_Unmelted_count = 0;
	long CELL_IS_POWDER_count = 0;
	long CELL_IS_LIQUID_count = 0;
	long CELL_IS_PARTIALLY_CRYSTALLIZED_count = 0;
	long UnmeltedPowder = 0;
	long CELL_IS_CRYSTALLIZED_count = 0;
	long CELL_IS_PARTICLE_count = 0;
	long CELL_IS_UNKNOWN_count = this->box.xPer* this->box.yPer* this->box.zPer;

	for (uint32_t z = z0; z < nz; z++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		//totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * RGB;
		totalOffset = (x0 + z * Nx * Ny + y0 * Nx);

		for (uint32_t y = y0; y < ny; y++)
		{
			long counterx = 0;
			MPI_Barrier(MPI_COMM_WORLD);

			for (uint32_t x = x0; x < nx; x++)
			{
				cellP currentCell = NULL;
				//unsigned char rgb[3] = { WHITE_RGB, WHITE_RGB, WHITE_RGB };

				//uint8_t state = CELL_IS_UNKNOWN;
				unsigned char state = CELL_IS_UNKNOWN;

				cellP* cp = getCellPerP(x, y, z);
				cellP c = *cp;
				//oriRepresentation oriDummy;

				if (!cellExists(c))
				{
					state = CELL_IS_POWDER;
					CELL_IS_POWDER_count++;
					CELL_IS_UNKNOWN_count--;

				}

				if (cellExists(c))
				{
					//cout << "exist" << endl;
					if (isCellRecrystallized(c))
					{
						state = CELL_IS_CRYSTALLIZED;
						CELL_IS_CRYSTALLIZED_count++;
						CELL_IS_UNKNOWN_count--;
					}

					else if (c->rxFrac == 0) {
						state = CELL_IS_POWDER;
						CELL_IS_POWDER_count++;
						CELL_IS_UNKNOWN_count--;
					}


					else if (InsufficientMelting(c)) {
						state = CELL_IS_UnmeltedPowder;
						CELL_IS_Unmelted_count++;
						CELL_IS_UNKNOWN_count--;
					}

					else
					{
						state = CELL_IS_PARTIALLY_CRYSTALLIZED;
						CELL_IS_PARTIALLY_CRYSTALLIZED_count++;
						CELL_IS_UNKNOWN_count--;
					}
				}

				xRowData[counterx] = state;
				counterx++;

			}
			//MPI_File_write_at(msFileHdl, totalOffset, xRowData, counterx * RGB, MPI_UNSIGNED_CHAR, &msFileStatus);
			//totalOffset += Nx * RGB;
			MPI_File_write_at(msFileHdl, totalOffset, xRowData, counterx, MPI_UNSIGNED_CHAR, &msFileStatus);
			totalOffset += Nx;

		}
	}

	delete[] CmsFileName;
	free(xRowData);

	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
	cout << "UnmeltedPowder" << "\t" << CELL_IS_Unmelted_count << endl;
	cout << "UNKNOWN" << "\t" << CELL_IS_UNKNOWN_count << endl;
	cout << "POWDER" << "\t" << CELL_IS_POWDER_count << endl;
	cout << "LIQUID" << "\t" << CELL_IS_LIQUID_count << endl;
	cout << "PARTIALLY_CRYSTALLIZED" << "\t" << CELL_IS_PARTIALLY_CRYSTALLIZED_count << endl;
	cout << "CRYSTALLIZED" << "\t" << CELL_IS_CRYSTALLIZED_count << endl;
	cout << "PARTICLE" << "\t" << CELL_IS_PARTICLE_count << endl;

	cout << " total_count" << CELL_IS_UNKNOWN_count + CELL_IS_POWDER_count + CELL_IS_LIQUID_count + CELL_IS_PARTIALLY_CRYSTALLIZED_count + CELL_IS_CRYSTALLIZED_count + CELL_IS_PARTICLE_count << endl;

}

// draw pole figure


void caHdl::draw_pole_figure(void)
{
	//TARGET FORMAT IS INT FOR PARAVIEW OR AVIZO ONLY RECRYSTALLIZED GRAINS (rxFrac < 1.0) ARE colored white
//CELL_IS_PARTICLE	0
//CELL_IS_INFECTED	1
//DEFORMED 			1+mydefgid
//RECRYSTALLIZED	1+mydefgid.size()+myrxgid
//<MPI_INT>uint

	//MPI_File msFileHdl;
	//MPI_Status msFileStatus;

	// create C-consistent file name for MPI I/O

	//stringstream msFileName;
	//msFileName << "SCORE."  "time." << this->time;
	//msFileName << "pole_figure.txt";

	string currentTime = std::to_string(this->time);
	string fileName = "SCORE." + currentTime + "pole_figure_z250.ctf";

	ofstream PoleFigureFile(fileName, ios::out);
	PoleFigureFile << "Channel Text File\n" <<
		"Prj	X30Mn22 ohne Rotation\n" <<
		"Author\n" <<
		"JobMode	Grid\n" <<
		"XCells	500\n" <<
		"YCells	500\n" <<
		"XStep	0.8\n" <<
		"YStep	0.8\n" <<
		"AcqE1	0.0000\n" <<
		"AcqE2	0.0000\n" <<
		"AcqE3	0.0000\n" <<
		"Euler angles refer to Sample Coordinate system(CS0)!Mag	80.0000	Coverage	99	Device	0	KV	20.0000	TiltAngle	70.0000	TiltAxis	0	DetectorOrientationE1	359.1714	DetectorOrientationE2	90.7337	DetectorOrientationE3	0.9590	WorkingDistance	17.9013	InsertionDistance	209.2969\n" <<
		"Phases	1\n" <<
		"2.866;2.866;2.866	90.000;90.000;90.000	Iron bcc(old)	11	229			J.Appl.Phys.[JAPIAU], Vol. 42, Seiten 4290 - 95\n"
		"3.660;3.660;3.660	90.000;90.000;90.000	Iron fcc	11	225			Z.Angew.Phys.[ZAPHAX], Vol. 23, Seiten 245 - 249\n"
		"5.108;6.777;4.540	90.000;90.000;90.000	Fe3C	3	62			J.Solid State Chem.[JSSCBI], Vol. 51, Seiten 246 - 252\n"
		"Phase	X	Y	Bands	Error	Euler1	Euler2	Euler3	MAD	BC	BS" << endl;

	
	
	//int msFileNameLength = msFileName.str().size();
	//char* CmsFileName = new char[msFileNameLength + 1];
	//strcpy(CmsFileName, msFileName.str().c_str());

	//working all nodes open the file in create and write-only mode
	//MPI_File_open(MPI_COMM_SELF, CmsFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl);

	//MK:: with 64-bit architectures, offsets calculations can be performed as longs 2^64 bit thus on petabyte files...
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized

	uint32_t x0 = this->x0;
	uint32_t y0 = this->y0;
	uint32_t z0 = this->z0;
	uint32_t nx = this->xPer + x0;
	uint32_t ny = this->yPer + y0;
	uint32_t nz = this->zPer + z0;
	uint32_t Nx = this->box.xPer;
	uint32_t Ny = this->box.xPer;

	//unsigned char* xRowData = NULL;
	//xRowData = new unsigned char[(nx - x0) * 12]; // 3 * (4 bytes) numbers 

	//long long totalOffset = 0;

	//MPI_File_seek(msFileHdl, totalOffset, MPI_SEEK_SET);

	//MPI_Barrier(MPI_COMM_WORLD);

	//unsigned char bytes[12] = { 0 };  // 3 numbers each with 4 bytes = 12 bytes

	for (uint32_t z = 250; z < 251; z++)
	{
		//MPI_Barrier(MPI_COMM_WORLD);
		//totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * 12;
		for (uint32_t y = 70; y < ny; y++)
		{
			//long counterx = 0;
			//MPI_Barrier(MPI_COMM_WORLD);
			for (uint32_t x = 0; x < nx; x++)
			{
				cellP currentCell = NULL;
				//uint32_t rgb[3] = { (uint32_t)0, (uint32_t)0, (uint32_t)0 };

				cellP* cp = getCellPerP(x, y, z);// z ihas been replaced with y
				cellP c = *cp;
				oriRepresentation oriDummy;

				if (!cellExists(c))
				{
					PoleFigureFile << "0\t" << fixed << setprecision(4) << z * 0.8000 << "\t" << y * 0.8000 << "\t" << "0\t3\t" << "0.0000" << "\t" << "0.0000" << "\t" << "0.0000" << "\t" << "0.0000\t145\t203" << endl;

				}

				if (cellExists(c))
				{
					if (isParticle(c))
					{

						//rgb[RED] = (uint32_t)BLACK;
						//rgb[GREEN] = (uint32_t)BLACK;
						//rgb[BLUE] = (uint32_t)BLACK;
					}
					else
					{
						grainP g = grainOfCell(c);
						long localOriIndex = g->oriIndex;
						oriRepresentationP localOri = grains.orientations[localOriIndex];
						oriDummy.copyOriFrom(localOri);

						//rgb[RED] = (uint32_t)(returnPositiveAngle(oriDummy.phi1) * 10000);
						//rgb[GREEN] = (uint32_t)(returnPositiveAngle(oriDummy.PHI) * 10000);
						//rgb[BLUE] = (uint32_t)(returnPositiveAngle(oriDummy.phi2) * 10000);
						PoleFigureFile << "1\t" << fixed << setprecision(4) << x * 0.8 << "\t" << y * 0.8 << "\t" << "8\t0\t" << oriDummy.phi1 * 180 / _PI_ << "\t" << oriDummy.PHI * 180 / _PI_ << "\t" << oriDummy.phi2 * 180 / _PI_ << "\t" << "0.4417\t145\t203" << endl;

					}
				}

				//for (int32_t i = 0; i < 3; i++)
				//	for (int32_t j = 3; j >= 0; j--)
				//	{
				//		//bytes[(3 - j) + 4 * i] = (rgb[i] >> (8 * j)) & 0xFF;// storing data as big endian form.
				//		bytes[j + 4 * i] = (rgb[i] >> (8 * j)) & 0xFF;// storing data as little endian form.
				//	}

				//memcpy(&xRowData[counterx * 12], bytes, 12 * sizeof(unsigned char));

				//counterx++;
			}

			/*			MPI_File_write_at(msFileHdl, totalOffset, xRowData, 12 * counterx, MPI_UNSIGNED_CHAR, &msFileStatus);
				totalOffset += ((long long)Nx * 12);*/
		}
	}





	// for yz plane
/*
	 
	string fileName1 = "SCORE.YZ" + currentTime + "pole_figure_z250.ctf";
	ofstream PoleFigureFile1(fileName, ios::out);
	PoleFigureFile1 << "Channel Text File\n" <<
		"Prj	X30Mn22 ohne Rotation\n" <<
		"Author\n" <<
		"JobMode	Grid\n" <<
		"XCells	500\n" <<
		"YCells	500\n" <<
		"XStep	0.8\n" <<
		"YStep	0.8\n" <<
		"AcqE1	0.0000\n" <<
		"AcqE2	0.0000\n" <<
		"AcqE3	0.0000\n" <<
		"Euler angles refer to Sample Coordinate system(CS0)!Mag	80.0000	Coverage	99	Device	0	KV	20.0000	TiltAngle	70.0000	TiltAxis	0	DetectorOrientationE1	359.1714	DetectorOrientationE2	90.7337	DetectorOrientationE3	0.9590	WorkingDistance	17.9013	InsertionDistance	209.2969\n" <<
		"Phases	1\n" <<
		"2.866;2.866;2.866	90.000;90.000;90.000	Iron bcc(old)	11	229			J.Appl.Phys.[JAPIAU], Vol. 42, Seiten 4290 - 95\n"
		"3.660;3.660;3.660	90.000;90.000;90.000	Iron fcc	11	225			Z.Angew.Phys.[ZAPHAX], Vol. 23, Seiten 245 - 249\n"
		"5.108;6.777;4.540	90.000;90.000;90.000	Fe3C	3	62			J.Solid State Chem.[JSSCBI], Vol. 51, Seiten 246 - 252\n"
		"Phase	X	Y	Bands	Error	Euler1	Euler2	Euler3	MAD	BC	BS" << endl;



	for (uint32_t z = 0; z < nz; z++)
	{
		//MPI_Barrier(MPI_COMM_WORLD);
		//totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * 12;
		for (uint32_t y = y0; y < ny; y++)
		{
			//long counterx = 0;
			//MPI_Barrier(MPI_COMM_WORLD);
			for (uint32_t x = 249; x < 250; x++)
			{
				cellP currentCell = NULL;
				//uint32_t rgb[3] = { (uint32_t)0, (uint32_t)0, (uint32_t)0 };

				cellP* cp = getCellPerP(x, y, z);// z ihas been replaced with y
				cellP c = *cp;
				oriRepresentation oriDummy;

				if (!cellExists(c))
				{
					PoleFigureFile1 << "0\t" << fixed << setprecision(4) << x * 0.8000 << "\t" << y * 0.8000 << "\t" << "0\t3\t" << "0.0000" << "\t" << "0.0000" << "\t" << "0.0000" << "\t" << "0.0000\t145\t203" << endl;

				}

				if (cellExists(c))
				{
					if (isParticle(c))
					{

						//rgb[RED] = (uint32_t)BLACK;
						//rgb[GREEN] = (uint32_t)BLACK;
						//rgb[BLUE] = (uint32_t)BLACK;
					}
					else
					{
						grainP g = grainOfCell(c);
						long localOriIndex = g->oriIndex;
						oriRepresentationP localOri = grains.orientations[localOriIndex];
						oriDummy.copyOriFrom(localOri);

						//rgb[RED] = (uint32_t)(returnPositiveAngle(oriDummy.phi1) * 10000);
						//rgb[GREEN] = (uint32_t)(returnPositiveAngle(oriDummy.PHI) * 10000);
						//rgb[BLUE] = (uint32_t)(returnPositiveAngle(oriDummy.phi2) * 10000);
						PoleFigureFile1 << "1\t" << fixed << setprecision(4) << x * 0.8 << "\t" << y * 0.8 << "\t" << "8\t0\t" << oriDummy.phi1 * 180 / _PI_ << "\t" << oriDummy.PHI * 180 / _PI_ << "\t" << oriDummy.phi2 * 180 / _PI_ << "\t" << "0.4417\t145\t203" << endl;

					}
				}

				//for (int32_t i = 0; i < 3; i++)
				//	for (int32_t j = 3; j >= 0; j--)
				//	{
				//		//bytes[(3 - j) + 4 * i] = (rgb[i] >> (8 * j)) & 0xFF;// storing data as big endian form.
				//		bytes[j + 4 * i] = (rgb[i] >> (8 * j)) & 0xFF;// storing data as little endian form.
				//	}

				//memcpy(&xRowData[counterx * 12], bytes, 12 * sizeof(unsigned char));

				//counterx++;
			}

			//		MPI_File_write_at(msFileHdl, totalOffset, xRowData, 12 * counterx, MPI_UNSIGNED_CHAR, &msFileStatus);
				//totalOffset += ((long long)Nx * 12);
		}
	}



	*/









	//delete[] CmsFileName;
	//delete[] xRowData;
	PoleFigureFile.close();
//	PoleFigureFile1.close();
	//MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}

void caHdl::write_prorosity(void)  //#as written for writng. declared in line 150 in caHdl.h 
{	
	//uint32_t x0 = this->x0;
	//uint32_t y0 = this->y0;
	//uint32_t z0 = this->z0;
	//uint32_t nx = this->xPer + x0;
	//uint32_t ny = this->yPer + y0;
	//uint32_t nz = this->zPer + z0;
	//uint32_t Nx = this->box.xPer;
	//uint32_t Ny = this->box.xPer;
	//Real count = 0;
	//
	// 
	ofstream  nomelt;
	nomelt.open("unmeltedposition.txt");
	nomelt << "slno \t\t\t tmax \t\t\t\ x \t\t\t y \t\t\t z \t\t\t" << endl;

	for (int i = 0; i < PT.tmax.size(); i++)


	{
		nomelt << i << " \t\t\t" << PT.tmax[i] << " \t\t\t" << PT.x[i]<<" \t\t\t"<< PT.y[i] << " \t\t\t"<< PT.z[i] << " \t\t\t"<<endl;

	}
	//for (uint32_t z = z0; z < nz; z++)
	//{
	//	//MPI_Barrier(MPI_COMM_WORLD);
	//	//totalOffset = (x0 + z * Nx * Ny + y0 * Nx) * 12;
	//	for (uint32_t y = y0; y < ny; y++)
	//	{
	//		//long counterx = 0;
	//		//MPI_Barrier(MPI_COMM_WORLD);
	//		for (uint32_t x = x0; x < nx; x++)
	//		{
	//			cellP currentCell = NULL;
	//			//uint32_t rgb[3] = { (uint32_t)0, (uint32_t)0, (uint32_t)0 };

	//			cellP* cp = getCellPerP(x, y, z);
	//			cellP c = *cp;
	//			if (cellExists(c))
	//			{
	//				if (InsufficientMelting(c))
	//				{
	//					//porosity k;
	//					long index = x * x + y * y +z *z;
	//					nomelt << x<< "\t\t\t" << y << "\t\t\t" << z <<"\t\t\t"<< endl; // nomelt[i].tmax << endl;
	//					count++;
	//					
	//				}
	//			}
	//		}
	//	}
	//}
	//nomelt << "total unmelted cells:" << count << endl;
	nomelt.close();
	


}
/*void caHdl::porosity(void)



*/
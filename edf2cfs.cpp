// A high-performance multi-threaded c++ CLI application to convert EDF (European Data Format) files to 
// Compressed Feature Set (CFS) format. The CFS format is used by the Z3Score sleep scoring system (https://z3score.com). 
// Instead of using polysomnography data in European Data Format (EDF, https://en.wikipedia.org/wiki/European_Data_Format), 
// the Z3Score system uses CFS files. CFS files are on an average 17X smaller than corresponding EDF files. 
// This reduces data overhead significantly. The format does not allow any user identifiable information ensuring anonymity.
//
// Patents pending (c)-2017 Amiya Patanaik amiyain@gmail.com
// This application has miltiple dependencies including 
// sigpack: C++ signal processing library http://sigpack.sourceforge.net/ 
// tclap: Templatized C++ Command Line Parser Library, 
// BOOST: C++ libraries 
// armadillo: C++ linear algebra library
// zlib:  DEFLATE compression algorithm
// Licensed under GPL v3

// edf2cfs.cpp : Defines the entry point for the console application.
//
extern "C" {
#include "edflib.h"
#include "zlib.h"
}
#include "sigpack/sigpack.h"
#include "tclap/CmdLine.h"
#include "tclap/ValueArg.h"
#include "SHA1.h"
#include "order32.h"
#include "resample.h"
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <bitset>
#include <cassert>
#include <sys/stat.h>
#include <string>
#include <thread>
#include <future> 
#include <functional>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <chrono>

#define FILTERORDER  (50)
#define SAMPLINGRATE  (100)
#define LITTLEENDIAN (O32_HOST_ORDER == O32_LITTLE_ENDIAN)
#define BOOST_FILESYSTEM_VERSION 3
#define BOOST_FILESYSTEM_NO_DEPRECATED 
#define DEBUG (false)
#define BR "<br />"

using namespace std;
using namespace arma;
namespace fs = ::boost::filesystem;

vec firBandPass(int N, double fl, double fh);
string removeExtension(const string& filename);
void writeByteReversed(FILE* file, Bytef* stream, int byteSize, unsigned long streamSize);
int roundInt(double r);
char *strlwr(char *str);

double findMultiplier(const string& units);

void showHeader(const char* filename, vector<string>& channelLabels);
bool convertFile(const char* filename, vector<string>* channelLabelsPtr, bool overwrite, ostringstream* streamMsgPointer);
void getAllFiles(const fs::path& root, const string& ext, vector<string>& filelist);

sp::FFTW fftObject(128);

int main(int argc, char *argv[]) {
	//Make sure IEEE-754 is supported
	assert(numeric_limits<float>::is_iec559 == true);

	//Number of Cores
	unsigned concurentThreadsSupported = thread::hardware_concurrency();

	if(concurentThreadsSupported==0)
		concurentThreadsSupported=2;

	vector<string> channelLabels;
	vector<string> filelist;
	string dirName;
	bool quiet;
	bool overwrite;
	bool saveLog;
	string logFile;
	ofstream lfile;

	//Parse command line

	try {
		TCLAP::CmdLine cmd("Usage: ./edf2cfs -a C3A2 -b C4A1 -x ELA2 -z ERA1 -q -o -l -d edfDir filename1.edf filename2.edf ... filenameN.edf\nIf no channels are given, then a selection menu will be shown.\n Use -d to provide a directory path with EDF files, -q to supress output, -o to overwrite and -l to save log.", ' ', "1.0");
		TCLAP::SwitchArg isquiet("q", "quiet", "silent mode", false);
		TCLAP::SwitchArg isoverwrite("o", "overwrite", "over write files", false);
		TCLAP::SwitchArg islog("l", "log", "save log", false);
		TCLAP::ValueArg<string> C3("a", "c3", "C3-A2 Channel Label", false, "NA", "C3-A2 Channel Label");
		TCLAP::ValueArg<string> C4("b", "c4", "C4-A1 Channel Label", false, "NA", "C4-A1 Channel Label");
		TCLAP::ValueArg<string> EL("x", "el", "EL-A2 Channel Label", false, "NA", "EL-A2 Channel Label");
		TCLAP::ValueArg<string> ER("z", "er", "ER-A1 Channel Label", false, "NA", "ER-A1 Channel Label");
		TCLAP::ValueArg<string> dir("d", "dir", "EDF Directory", false, "NA", "EDF Directory");
		TCLAP::UnlabeledMultiArg<string> files("filenames", "List of filename", false, "List of EDF files", false);

		cmd.add(files);
		cmd.add(C3);
		cmd.add(C4);
		cmd.add(EL);
		cmd.add(ER);
		cmd.add(dir);
		cmd.add(isquiet);
		cmd.add(isoverwrite);
		cmd.add(islog);

		if (argc < 2) {
			cout << "No EDF files provided\n";
			cout << "./edf2cfs -h for usage details.\n";
			std::system("read -n 1 -s -p \"Press any key to continue...\"");
			return(1);
		}
		
		cmd.parse(argc, argv);
		channelLabels.push_back(C3.getValue());
		channelLabels.push_back(C4.getValue());
		channelLabels.push_back(EL.getValue());
		channelLabels.push_back(ER.getValue());
		dirName = dir.getValue();
		quiet = isquiet.getValue();
		overwrite = isoverwrite.getValue();
		saveLog = islog.getValue();
		filelist = files.getValue();
		if(strcmp(dirName.c_str(),"NA") != 0){
			fs::path dirPath{dirName};
			getAllFiles(dirPath,".edf",filelist);
		}

		if(filelist.empty()){
			cout << "No EDF files found.\n";
			cout << "./edf2cfs -h for usage details.\n";
			std::system("read -n 1 -s -p \"Press any key to continue...\"");
			return(1);
		}
	}
	catch (TCLAP::ArgException& e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		return(1);
	}

	//Check if channel numbers are provided
	if (strcmp(channelLabels[0].c_str(),"NA") == 0 || strcmp(channelLabels[1].c_str(), "NA") == 0  || strcmp(channelLabels[2].c_str(), "NA") == 0  || strcmp(channelLabels[3].c_str(), "NA") == 0 ) {
		//Channels not provided 
		showHeader(filelist[0].c_str(), channelLabels);
	}

	//if logging is on
	if(saveLog){
		const std::locale fmt(std::locale::classic(), new boost::posix_time::time_facet("%d-%b-%Y-%H%M"));
		boost::posix_time::ptime tnow(boost::posix_time::second_clock::local_time());
		std::ostringstream os;
		os.imbue(fmt);
		os << tnow;
		fs::path filePath{filelist[0]};
		fs::path fullPath = fs::complete(filePath);
		fs::path basePath = fullPath.parent_path();
		logFile = basePath.string() + "/" + os.str() + "_log.html";
		cout<<"Log will be saved at:\n" << logFile << endl;
		try{
			lfile.open(logFile,ios::out);
		}
		catch (ifstream::failure e) {
			saveLog = false;
			lfile.close();
		}



		if(saveLog) {
			lfile << "<!doctype html>\n<html lang='en'>\n<head>\n" 
			"<meta charset='utf-8'>\n\n  <title>EDF to CFS Log</title>\n"  
			"<meta name='description' content='Conversion Log'>\n" 
			"<meta name='author' content='Amiya Patanaik'>\n\n  </head>\n\n<body>\n";

			lfile << "<p>Logging Started at: " << os.str() << BR <<endl;
			lfile << "C3-A2 Channel Label: " << channelLabels[0] << BR <<endl;
			lfile << "C4-A1 Channel Label: " << channelLabels[1] << BR <<endl;
			lfile << "EL-A2 Channel Label: " << channelLabels[2] << BR <<endl;
			lfile << "ER-A1 Channel Label: " << channelLabels[3] << BR <<endl;
			lfile << "</p><hr>" <<endl; 
		}

	}

	//Start conversion 
	vec successCounter = zeros<vec>(filelist.size());
	printf("Processing upto %d files simultanously...\n", concurentThreadsSupported);
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < (int)filelist.size(); i += concurentThreadsSupported) {

		vector<future<bool>> workers;
		vector<ostringstream> streamMsg(concurentThreadsSupported);

		workers.reserve(concurentThreadsSupported);

		for (int j = 0; j< concurentThreadsSupported; j++) {

			if (i + j >= (int)filelist.size())
				break;

			workers.push_back(async(convertFile, filelist[i+j].c_str(), &channelLabels, overwrite, &streamMsg[j]));

		}

		//Now wait
		for (int j = 0; j < workers.size(); j++) {
			successCounter(i + j) = workers[j].get();
			if (successCounter(i + j) == 0) //If failed always print output
				//cout << streamMsg[j].str();
				if(!saveLog){
					cout << "ERROR: Filename: " << filelist[i+j].c_str() << ", please enable logging to see details.\n";
				}
				else {
					cout << "ERROR: Filename: " << filelist[i+j].c_str() << ", please check log.\n";
				}
				
			else if (!quiet)
				//cout << streamMsg[j].str();
				cout << "Filename: " << filelist[i+j].c_str() << ", processed successfully\n";

			if(saveLog){
				lfile << streamMsg[j].str();
			}
		}
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto intms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	int intSecs = (int)(intms.count()/1000);

	printf("%lu Files processed in %d seconds.\n%d Files converted successfully. %lu Files could not be converted.\n", filelist.size(),intSecs,(int)accu(successCounter), filelist.size()- (int)accu(successCounter));

	if(saveLog){
		lfile << filelist.size() << " Files processed in " << intSecs<< " seconds." << BR <<endl;
		lfile << (int)accu(successCounter) << " Files converted successfully. " << (filelist.size()- (int)accu(successCounter)) << " Files could not be converted.<br />";
		lfile.close();
	}
	std::system("read -n 1 -s -p \"Press any key to continue...\"");
}

void showHeader(const char* filename, vector<string>& channelLabels) {

	struct edf_hdr_struct hdr;
	int nC3, nC4, nEL, nER;

	if (edfopen_file_readonly(filename, &hdr, EDFLIB_READ_ALL_ANNOTATIONS)) {
		switch (hdr.filetype) {
		case EDFLIB_MALLOC_ERROR: printf("\nMemory Error.\n\n");
			break;
		case EDFLIB_NO_SUCH_FILE_OR_DIRECTORY: printf("\nCan not open file, no such file or directory\n\n");
			break;
		case EDFLIB_FILE_CONTAINS_FORMAT_ERRORS: printf("\nThe file is not EDF(+) or BDF(+) compliant\n"
			"(it contains format errors)\n\n");
			break;
		case EDFLIB_MAXFILES_REACHED: printf("\nToo many files opened\n\n");
			break;
		case EDFLIB_FILE_READ_ERROR: printf("\nA read error occurred\n\n");
			break;
		case EDFLIB_FILE_ALREADY_OPENED: printf("\nFile has already been opened\n\n");
			break;
		default: printf("\nUnknown error\n\n");
			break;
		}

		exit(EXIT_FAILURE);
	}

	int Nmax = hdr.edfsignals;
	int hdl = hdr.handle;

	printf("Please make sure all files share the same channel labels.\n");
	printf("Following channels are found:\n");

	for (int i = 0; i<hdr.edfsignals; i++) {
		printf("%d: %s\n", i + 1, hdr.signalparam[i].label);
	}

	printf("Please select the C3:A2 channel number: \n");
	cin >> nC3;
	printf("Please select the C4:A1 channel number: \n");
	cin >> nC4;
	printf("Please select the EOGl:A2 channel number: \n");
	cin >> nEL;
	printf("Please select the EOGr:A1 channel number: \n");
	cin >> nER;

	nC3 -= 1;
	nC4 -= 1;
	nEL -= 1;
	nER -= 1;

	if (nC3 > Nmax || nC4 > Nmax || nEL > Nmax || nER > Nmax || nC3 < 0 || nC4 < 0 || nEL < 0 || nER < 0) {
		printf("Invalid Channel Number.\n");
		edfclose_file(hdl);
		exit(EXIT_FAILURE);
	}

	channelLabels[0] = (strlwr(hdr.signalparam[nC3].label));
	channelLabels[1] = (strlwr(hdr.signalparam[nC4].label));
	channelLabels[2] = (strlwr(hdr.signalparam[nEL].label));
	channelLabels[3] = (strlwr(hdr.signalparam[nER].label));

	edfclose_file(hdl);

}

bool convertFile(const char* filename, vector<string>* channelLabelsPtr, bool overwrite, ostringstream* streamMsgPointer) {

	ostringstream &streamMsg = *(streamMsgPointer);
	vector<string> &channelLabels = *(channelLabelsPtr);
	streamMsg << "<p>Filename: " << filename << BR <<endl;

	//filename for CFS file.
	string baseName = removeExtension(string(filename)) + ".cfs";

	if (!overwrite) {
		if (fs::exists(baseName)) {

			streamMsg << "<strong style='color:red;'>ERROR: File already converted.</strong><br /></p>\n";
			return false;
		}
	}
	struct edf_hdr_struct hdr;
	double fC3, fC4, fEL, fER;
	int nC3, nC4, nEL, nER;
	long long samplesRead;

	vector<string> allLabels;

	if (edfopen_file_readonly(filename, &hdr, EDFLIB_READ_ALL_ANNOTATIONS)) {

		switch (hdr.filetype) {
		case EDFLIB_MALLOC_ERROR: streamMsg << "<strong style='color:red;'>ERROR: Memory Error.</strong><br />\n\n</p>";
			break;
		case EDFLIB_NO_SUCH_FILE_OR_DIRECTORY: streamMsg << "<strong style='color:red;'>ERROR: Can not open file, no such file or directory</strong><br />\n\n</p>";
			break;
		case EDFLIB_FILE_CONTAINS_FORMAT_ERRORS: streamMsg << "<strong style='color:red;'>ERROR: The file is not EDF(+) or BDF(+) "
		"compliant (it contains format errors)</strong><br />\n\n</p>";
			break;
		case EDFLIB_MAXFILES_REACHED: streamMsg << "<strong style='color:red;'>ERROR: Too many files opened</strong><br />\n\n</p>";
			break;
		case EDFLIB_FILE_READ_ERROR: streamMsg << "<strong style='color:red;'>ERROR: A read error occurred</strong><br />\n\n</p>";
			break;
		case EDFLIB_FILE_ALREADY_OPENED: streamMsg << "<strong style='color:red;'>ERROR: File has already been opened</strong><br />\n\n</p>";
			break;
		default: streamMsg << "<strong style='color:red;'>ERROR: Unknown error</strong><br />\n\n</p>";
			break;
		}

		return false;
	}

	int Nmax = hdr.edfsignals;
	int hdl = hdr.handle;

	for (int i = 0; i<hdr.edfsignals; i++) {
		allLabels.push_back(strlwr(hdr.signalparam[i].label));
	}

	vector<string>::iterator it;

	//C3 channel
	it = find(allLabels.begin(), allLabels.end(), channelLabels[0]);
	if (it != allLabels.end()) {
		nC3 = distance(allLabels.begin(), it);
	}
	else {

		streamMsg << "<strong style='color:red;'>Error: C3 label not found!</strong><br />\n</p>" << endl;
		edfclose_file(hdl);
		return false;
	}

	//C4 channel
	it = find(allLabels.begin(), allLabels.end(), channelLabels[1]);
	if (it != allLabels.end()) {
		nC4 = distance(allLabels.begin(), it);
	}
	else {

		streamMsg << "<strong style='color:red;'>Error: C4 label not found!</strong><br />\n</p>" << endl;
		edfclose_file(hdl);
		return false;
	}

	//EL channel
	it = find(allLabels.begin(), allLabels.end(), channelLabels[2]);
	if (it != allLabels.end()) {
		nEL = distance(allLabels.begin(), it);
	}
	else {

		streamMsg << "<strong style='color:red;'>Error: EL label not found!</strong><br />\n</p>" << endl;
		edfclose_file(hdl);
		return false;
	}

	//ER channel
	it = find(allLabels.begin(), allLabels.end(), channelLabels[3]);
	if (it != allLabels.end()) {
		nER = distance(allLabels.begin(), it);
	}
	else {

		streamMsg << "<strong style='color:red;'>Error: ER label not found!</strong><br />\n</p>" << endl;
		edfclose_file(hdl);
		return false;
	}

	//Read sampling rates 
	fC3 = ((double)hdr.signalparam[nC3].smp_in_datarecord / (double)hdr.datarecord_duration) * EDFLIB_TIME_DIMENSION;
	fC4 = ((double)hdr.signalparam[nC4].smp_in_datarecord / (double)hdr.datarecord_duration) * EDFLIB_TIME_DIMENSION;
	fEL = ((double)hdr.signalparam[nEL].smp_in_datarecord / (double)hdr.datarecord_duration) * EDFLIB_TIME_DIMENSION;
	fER = ((double)hdr.signalparam[nER].smp_in_datarecord / (double)hdr.datarecord_duration) * EDFLIB_TIME_DIMENSION;

	//Read measurement units
	string fC3unit = hdr.signalparam[nC3].physdimension;
	string fC4unit = hdr.signalparam[nC4].physdimension;
	string fELunit = hdr.signalparam[nEL].physdimension;
	string fERunit = hdr.signalparam[nER].physdimension;

	//Ensure units are in uV
	double fC3mult = findMultiplier(fC3unit);
	double fC4mult = findMultiplier(fC4unit);
	double fELmult = findMultiplier(fELunit);
	double fERmult = findMultiplier(fERunit);

	if (fC3mult < 0 || fC4mult < 0 || fELmult < 0 || fERmult < 0) {

		streamMsg << "<strong style='color:red;'>ERROR: Invalid measurement unit. (must be nV, uV, mV or V)</strong><br />\n</p>" << endl;
		edfclose_file(hdl);
		return false;
	}


	if ((int)fC3 != (int)fC4) {

		streamMsg << "<strong style='color:red;'>Error: C3 and C4 sampling rates must be same.</strong><br />\n</p>" << endl;
		edfclose_file(hdl);
		return false;
	}

	long long totSamples = hdr.signalparam[nC3].smp_in_file;


	vector<double> bufC3(totSamples);
	vector<double> bufC4(totSamples);
	vector<double> bufEL(totSamples);
	vector<double> bufER(totSamples);

	streamMsg << "Total Samples found: " << totSamples << BR << endl;
	streamMsg << "C3:A2 channel, sampling rate: " << fC3 << "Hz measured in " << fC3unit << BR << endl;


	edfrewind(hdl, nC3);
	samplesRead = edfread_physical_samples(hdl, nC3, totSamples, &bufC3[0]);

	if (samplesRead == (-1)) {

		streamMsg << "\n<strong style='color:red;'>ERROR: reading channel C3 data.</strong><br />\n</p>\n";
		edfclose_file(hdl);
		return false;
	}

	streamMsg << "C4:A1 channel, sampling rate: " << fC4 << "Hz measured in " << fC4unit << BR << endl;

	edfrewind(hdl, nC4);
	samplesRead = edfread_physical_samples(hdl, nC4, totSamples, &bufC4[0]);

	if (samplesRead == (-1)) {

		streamMsg << "\n<strong style='color:red;'>ERROR: reading channel C4 data.</strong><br />\n</p>\n";
		edfclose_file(hdl);
		return false;
	}

	streamMsg << "EOGl:A2 channel, sampling rate: " << fEL << "Hz measured in " << fELunit << BR << endl;

	edfrewind(hdl, nEL);
	samplesRead = edfread_physical_samples(hdl, nEL, totSamples, &bufEL[0]);

	if (samplesRead == (-1)) {
		streamMsg << "\n<strong style='color:red;'>ERROR: reading channel EOG-l  data.</strong><br />\n</p>\n";
		edfclose_file(hdl);
		return false;
	}


	streamMsg << "EOGr:A1 channel, sampling rate: " << fER << "Hz measured in " << fERunit << BR << endl;

	edfrewind(hdl, nER);
	samplesRead = edfread_physical_samples(hdl, nER, totSamples, &bufER[0]);

	if (samplesRead == (-1)) {

		streamMsg << "\n<strong style='color:red;'>ERROR: reading channel EOG-R data.</strong><br />\n</p>\n";
		edfclose_file(hdl);
		return false;
	}


	edfclose_file(hdl);

	//Convert to Armadillo vectors
	//This is in-place conversion to avoid copy
	//Please remember!!
	vec dataC3(&bufC3[0], totSamples, false, true);
	vec dataC4(&bufC4[0], totSamples, false, true);
	vec dataEL(&bufEL[0], totSamples, false, true);
	vec dataER(&bufER[0], totSamples, false, true);

	if(DEBUG){
		dataC3.save("C3_Orig.csv", arma_ascii);
	}


	//Initialize Order 50 FIR bandpass filter weights
	vec filterEEG = firBandPass(FILTERORDER, 0.3 * 2 / fC3, 45 * 2 / fC3);
	vec filterEOGL = firBandPass(FILTERORDER, 0.3 * 2 / fEL, 12 * 2 / fEL);
	vec filterEOGR;

	if (fER == fEL)
		filterEOGR = firBandPass(FILTERORDER, 0.3 * 2 / fER, 12 * 2 / fER);
	else
		filterEOGR = filterEOGL;

	//mean computation and FIR filtering
	//vec eeg = conv((dataC3*fC3mult + dataC4*fC4mult) / 2.0, filterEEG, "same");
	vec eeg = (conv(dataC3*fC3mult, filterEEG, "same") + conv(dataC4*fC4mult, filterEEG, "same"))/2.0;
	vec eogl = conv(dataEL*fELmult, filterEOGL, "same");
	vec eogr = conv(dataER*fERmult, filterEOGR, "same");

	if(DEBUG){
		filterEEG.save("EEGb.csv", arma_ascii);
		eeg.save("eegFiltered.csv", arma_ascii);
	}


	vec eegFilt;
	vec eoglFilt;
	vec eogrFilt;

	//downsampling eeg signal to 100Hz
	if ((int)fC3 != SAMPLINGRATE) {
		vector<double> input = arma::conv_to< vector<double> >::from(eeg);
		vector<double> output;
		resample(SAMPLINGRATE,fC3,input,output);
		eegFilt = arma::conv_to<vec>::from(output);
	}
	else {
		eegFilt = eeg;
	}

	//downsampling eog signals to 100Hz
	if ((int)fEL != SAMPLINGRATE) {
		vector<double> input = arma::conv_to< vector<double> >::from(eogl);
		vector<double> output;
		resample(SAMPLINGRATE,fEL,input,output);
		eoglFilt = arma::conv_to<vec>::from(output);
	}
	else {
		eoglFilt = eogl;
	}

	if ((int)fER != SAMPLINGRATE) {
		vector<double> input = arma::conv_to< vector<double> >::from(eogr);
		vector<double> output;
		resample(SAMPLINGRATE,fER,input,output);
		eogrFilt = arma::conv_to<vec>::from(output);
	}
	else {
		eogrFilt = eogr;
	}

	if(DEBUG){
		eegFilt.save("EEGresampled.csv", arma_ascii);
		eoglFilt.save("EOGlresampled.csv", arma_ascii);
		eogrFilt.save("EOGrresampled.csv", arma_ascii);
	}


	//compute spectrogram
	int epochs = (int)(((double)eegFilt.n_elem) / 3000);

	static int epochSize = 32 * 32 * 3;
	vec payload(epochs*epochSize, fill::zeros);

	vec hamWindow = sp::hamming(128);

	for (int i = 0; i<epochs; ++i) {
		for (int j = 0; j<3000 - 128; j += 90) {

			//Time index
			int tIDX = (int)j / 90;

			vec eegEpochData = eegFilt.subvec(i * 3000 + j, i * 3000 + j + 128 - 1)%hamWindow;
			vec eegFFTResult = abs(fftObject.fft(eegEpochData));

			vec eoglEpochData = eoglFilt.subvec(i * 3000 + j, i * 3000 + j + 128 - 1)%hamWindow;
			vec eoglFFTResult = abs(fftObject.fft(eoglEpochData));

			vec eogrEpochData = eogrFilt.subvec(i * 3000 + j, i * 3000 + j + 128 - 1)%hamWindow;
			vec eogrFFTResult = abs(fftObject.fft(eogrEpochData));

			payload.subvec(i*epochSize + tIDX * 32, i*epochSize + tIDX * 32 + 32 - 1) = eegFFTResult.subvec(0, 32 - 1);
			payload.subvec(i*epochSize + 32 * 32 + tIDX * 32, i*epochSize + 32 * 32 + tIDX * 32 + 32 - 1) = eoglFFTResult.subvec(0, 32 - 1);
			payload.subvec(i*epochSize + 32 * 32 * 2 + tIDX * 32, i*epochSize + 32 * 32 * 2 + tIDX * 32 + 32 - 1) = eogrFFTResult.subvec(0, 32 - 1);


		}
	}

	if(DEBUG){
		payload.save("payload.csv", arma_ascii);
	}


	//Convert Double to Float (IEEE-754) to save space
	vector<float> fPayload = conv_to< vector<float> >::from(payload);

	//Convert float stream to binary stream
	unsigned long sourceLen = (fPayload.size())*sizeof(float);
	Bytef* istream = reinterpret_cast<Bytef*>(&fPayload[0]);

	//HEADER
	char signature[] = { 'C','F','S' };
	uint8_t version = 1;
	uint8_t nFreq = 32;
	uint8_t nTimes = 32;
	//next 16 bits
	uint8_t nChannels = 3;
	uint16_t nEpochs = epochs;
	uint8_t compression = true;
	uint8_t hash = true;


	//Find SHA1 hash of stream
	CSHA1 sha1;
	sha1.Update(istream, sourceLen);
	sha1.Final();
	Bytef* SHAdigest = new Bytef[20];
	if (!sha1.GetHash(SHAdigest)) {

		streamMsg << "<strong style='color:red;'>ERROR: Problem in conversion! SHA1 Failed...</strong><br />\n</p>\n";
		return false;
	}

	unsigned long destLen = compressBound(sourceLen);

	Bytef* ostream = new Bytef[destLen];

	int res = compress(ostream, &destLen, istream, sourceLen);

	if (res == Z_BUF_ERROR) {

		streamMsg << "<strong style='color:red;'>ERROR: Buffer was too small!</strong><br />\n</p>\n";
		return false;
	}
	if (res == Z_MEM_ERROR) {

		streamMsg << "<strong style='color:red;'>ERROR: Not enough memory for compression!</strong><br />\n</p>\n";
		return false;
	}

	FILE* file = fopen(baseName.c_str(), "wb");
	if (!file) {

		streamMsg << "<strong style='color:red;'>ERROR: Opening" << filename << "</strong><br />\n</p>" <<endl;
		return false;
	}
	//Start writing 
	//Header 11 bytes
	fwrite((Bytef*)signature, 1, 3, file);
	fwrite((Bytef*)&version, 1, 1, file);
	fwrite((Bytef*)&nFreq, 1, 1, file);
	fwrite((Bytef*)&nTimes, 1, 1, file);
	fwrite((Bytef*)&nChannels, 1, 1, file);
	if (LITTLEENDIAN)
		fwrite((Bytef*)&nEpochs, 1, 2, file);
	else
		writeByteReversed(file, (Bytef*)&nEpochs, 2, 2);
	fwrite((Bytef*)&compression, 1, 1, file);
	fwrite((Bytef*)&hash, 1, 1, file);
	//SHA1 20 bytes
	if (LITTLEENDIAN)
		fwrite(SHAdigest, 1, 20, file);
	else
		writeByteReversed(file, SHAdigest, 20, 20);
	//stream
	if (LITTLEENDIAN)
		fwrite(ostream, 1, destLen, file);
	else
		writeByteReversed(file, ostream, 4, destLen);
	//done
	fclose(file);

	streamMsg << "\n</p>";


	return(true);

}




vec firBandPass(int N, double fl, double fh) {
	vec b(N+1), h(N+1);
    h = sp::hamming(N+1);
    for (int i=0;i<N+1;i++)   {
        b[i] = h[i]*( sp::sinc(fh*(i-N/2.0))*fh - sp::sinc(fl*(i-N/2.0))*fl );
    }
    return b;
}

string removeExtension(const string& filename) {
	size_t lastdot = filename.find_last_of(".");
	if (lastdot == string::npos) return filename;
	return filename.substr(0, lastdot);
}

void writeByteReversed(FILE* file, Bytef* stream, int byteSize, unsigned long streamSize) {
	for (unsigned long offset = 0; offset < streamSize; offset += byteSize) {
		for (int j = byteSize - 1; j < 0; j--) {
			std::fwrite(stream + offset + j, 1, 1, file);
		}
	}
}

int roundInt(double r) {
	return (r > 0.0) ? (r + 0.5) : (r - 0.5);
}

char *strlwr(char *str) {
	unsigned char *p = (unsigned char *)str;

	while (*p) {
		*p = tolower(*p);
		p++;
	}

	return str;
}


double findMultiplier(const string& units) {

	if (units.compare(0,2,"nV") == 0)
		return 0.001;
	else if (units.compare(0, 2, "uV") == 0)
		return 1.0;
	else if (units.compare(0, 2, "mV") == 0)
		return 1000;
	else if (units.compare(0, 1, "V") == 0)
		return 1000000;
	else
		return -1;

}

void getAllFiles(const fs::path& root, const string& ext, vector<string>& filelist){
	if(!fs::exists(root) || !fs::is_directory(root)) return;

    fs::directory_iterator it(root);
    fs::directory_iterator endit;

    while(it != endit)
    {
        if(fs::is_regular_file(*it) && it->path().extension() == ext) filelist.push_back(it->path().string());
        ++it;

    }


}
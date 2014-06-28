// Copyright (c) 2000-2008, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are not permitted without specific prior written permission
// from TEES.

// If written permission is obtained for redistribution or further use, the
// following conditions must be met:

// 1) Redistributions of source code must retain the above copyright notice,
// this list of conditions and the disclaimer below.

// 2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions, and the disclaimer below in the documentation and/or
// other materials provided with the distribution.

// 3) Neither the name of TEES, the name of the Texas A&M University System, nor
// the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
#include "Tricall_Utility.h"
#define COMMENT_CHAR '#'

using namespace std;
using std::cout;
using std::ofstream;
using std::ios;

////////////////////////////////////////////////////////////////////////////////
///
/// @brief .
///
/// This class is used 
////////////////////////////////////////////////////////////////////////////////
class Tricall_Pre_Mesh_Gen
{
private:
	//

public:
	Tricall_Pre_Mesh_Gen(void);
	~Tricall_Pre_Mesh_Gen(void);

	void getlinesc(std::ifstream& inputfile, std::string& line);
	double String2Double(std::string s);
	void Tokenize(std::string& line, std::vector<std::string>& tokens,
		std::string& delimiters);
	double Xintersect(double x1, double y1, double x2, double y2, double y3);
	double Yintersect(double x1, double y1, double x2, double y2, double x3);
	template <class T> void printVector(vector <T> v);
	template <class T> void printVectorOfVector(vector < vector < T > > v);
	void trim(string&s);
	void SegmentSplitX(vector <double> & X, vector< vector <double> > & seg2point, vector < double > & xcoord,
		vector <double> & ycoord, vector< vector <double> > & ptatt, vector<double> & ptbdm, int & numseg, 
		int & numpts);
	void SegmentSplitY(vector <double> & Y, vector< vector <double> > & seg2point, vector < double > & xcoord, vector <double> & ycoord, vector< vector <double> > & ptatt, vector<double> & ptbdm, int & numseg, int & numpts);	
	void Grid(vector <double> & X, vector <double> & Y, vector <double> & xcoord, vector <double> & ycoord,
		vector < vector <double> > & ptatt, vector <double> & ptbdm, vector < vector <double> > & seg2point, 
		int & numatt, int & numpts, int & numbdm, int & numseg, int & numbdmseg);
	vector < vector < vector <double> > > PointSplit(vector <double> xcoord, vector <double> ycoord,
		vector <double> X, vector <double> Y, vector <vector <double> > ptatt, vector <double> ptbdm, 
		int numpts, int numatt, int numbdm);	
	vector < vector < vector <double> > > SegSplit(vector < vector < vector <double> > > polycell, 
		vector < vector <double> > seg2point, vector <double> X, vector <double> Y, int numbdm, int numseg);
	void Renumber (vector <double> & X, vector <double> & Y, vector < vector < vector <double> > > & polycell,
		vector < vector < vector <double> > > & segcell);
	int getpolyfiles(int, char**);
	int getpolyfiles_back();
};


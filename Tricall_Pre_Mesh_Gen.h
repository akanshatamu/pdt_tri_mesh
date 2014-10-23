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
#include <exception>
#include <stapl/runtime/runtime.hpp>
#include <stapl/runtime/executor/anonymous_executor.hpp>
#include <stapl/runtime/rmi/async_rmi.hpp>
#include <stapl/runtime/p_object.hpp>
#include <stapl/runtime/rmi_handle.hpp>
#include <iterator>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>
#include "Tricall_Utility.h"
#define COMMENT_CHAR '#'

using namespace std;
using std::cout;
using std::ofstream;
using std::ios;
using std::map;

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
  void Grid(vector <double> & X, vector <double> & Y, vector <double> & xcoord, vector <double> & ycoord,vector <vector<double> > & ptatt, vector <double> & ptbdm, vector <vector <double> > & seg2point, int & numatt, int & numbdm,int & numbdmseg);
  vector < vector < vector <double> > > PointSplit(vector <double> xcoord, vector <double> ycoord,
    vector <double> X, vector <double> Y, vector <vector <double> > ptatt, vector <double> ptbdm, 
    int numpts, int numatt, int numbdm);	
  vector < vector < vector <double> > > SegSplit(vector < vector < vector <double> > > polycell, 
    vector < vector <double> > seg2point, vector <double> X, vector <double> Y, int numbdm, int numseg);
  void Renumber (vector <double> & X, vector <double> & Y, vector < vector < vector <double> > > & polycell,
    vector < vector < vector <double> > > & segcell);
  int generateMesh(); 
  void getBoundaryNodes(vector <double>& , int bdrIndex, double bdrVal, vector <double>&,int,int,int,int,int,vector <int>&);
  void updateContainers(map<int,double>&,vector <double>&,vector <int>&);
  void load_geom(int & ntri, vector <double> & AttID, vector < vector <int> > & tri2edge,vector < vector <int> > & edge2tri, vector < vector <int> > & edge2pt, vector < vector <int> > & tri2pt,vector < vector <double> > & pt2coord, string & fileprefix, int & index);
  vector <vector<double> > SortSegments(vector <vector<double> > potentials, vector <double> ycoord);
  int WhichSide(double Xmin, double Xmax, double Ymin, double Ymax, double x, double y);
  vector <int> WhichTwo(double Xmin,double Xmax,double Ymin,double Ymax,double x1,double x2,double y1,double y2);
  int IsInside(double x, double y, double Xmin, double Xmax, double Ymin, double Ymax);
  void SegmentSplit(string & filename, double & Xmin, double & Xmax, double & Ymin, double & Ymax, vector <vector <double> > & seg2point, vector <double> & xcoord, vector <double> & ycoord, vector <vector <double> > & ptatt, vector <double> & ptbdm, int & numatt, int & numbdm, int & numbdmseg, vector <double> & newxcoord, vector <double> & newycoord, vector <vector <double> > & newptatt, vector <double> & newptbdm, vector <vector <double> > & SubDomSegments);
  int Neigh(int NeighX, int NeighY, map <string, int> SubCoordToInt);
  // Sorts based on the second column in a vector of vectors
    static bool sortFunc(const vector <double> &p1, const vector <double> &p2)
    {return p1[1] <= p2[1];}
    void Centroid(vector <vector<double> > & polygonpts, double & Cx, double & Cy);
  vector <vector <double> > RegionalAtt(string fileprefix,double Xmin,double Xmax,double Ymin,double Ymax);
}; 


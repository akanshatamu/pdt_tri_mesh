#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "tricall_utility.h"

using namespace std;

#define COMMENT_CHAR '#'

//utility functions
 inline void getlinesc(std::ifstream& inputfile, std::string& line)
 {
	std::getline(inputfile, line);
	if(line[0] == COMMENT_CHAR)
	getlinesc(inputfile, line);
 }
inline double String2Double(std::string s)
{
  double val;
  std::istringstream iss(s);
  iss >> val;
  return val;
}
inline void Tokenize(std::string& line, std::vector<std::string>& tokens,
  std::string& delimiters)
{
  std::string::size_type lastPos = line.find_first_not_of(delimiters, 0);
  std::string::size_type pos = line.find_first_of(delimiters, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    tokens.push_back(line.substr(lastPos, pos - lastPos));
    lastPos = line.find_first_not_of(delimiters, pos);
    pos = line.find_first_of(delimiters, lastPos);
  }
}

// Given two points that form a line (input as a 1x2 [x y] vector) and a y-plane that intersects 
// the segment formed by the two points, Xintersect finds the x-coordinate where the segment and 
// plane intersect
double Xintersect(double x1, double y1, double x2, double y2, double y3)
{
	
	// slope of the segment
	double m = (y2-y1)/(x2-x1);

	double x3 = ((y3-y1)/m) + x1;

	return x3;
}

// Given that two points form a line, and a x-plane (x3) intersects that segment,
// this function finds the y coordinate where the line and plane intersect
double Yintersect(double x1, double y1, double x2, double y2, double x3)
{
	
	//slope of the segment
	double m = (y2-y1)/(x2-x1);

	double y3 = (m*(x3-x1)) + y1;

	return y3;
}

template <class T> void printVector(vector <T> v)
{
	for (int i = 0; i < v.size(); i = i + 1)
	{
		cout << "Element " << i << " = " << v[i] << endl;
	}
}
template <class T> void printVectorOfVector(vector < vector < T > > v)
{
	for (int i = 0; i < v.size(); i = i + 1)
	{
		cout << "Element "<< i << " = ";
		for (int j = 0; j < v[i].size(); j = j + 1)
		{
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
}
void trim(string&s)
{
	size_t p = s.find_first_not_of("\t");
	s.erase(0,p);

	p = s.find_last_not_of("\t");
	if(string::npos != p)
		s.erase(p+1);
}

//Added by AKANSHA
//Call triangle library to generate triangle meshes
//Input: Poly file
//Output:- Triangle mesh related files
void callTriangleLibrary(string &s_poly_name){

}

// splitting the segments according the user-defined X sections
void SegmentSplitX(vector <double> & X, vector< vector <double> > & seg2point, vector < double > & xcoord,
	vector <double> & ycoord, vector< vector <double> > & ptatt, vector<double> & ptbdm, int & numseg, 
	int & numpts)
{
	// size of the X vector
	int Xlength = X.size();
	// number of attributes
	int numatt = ptatt[1].size();
	// number of boundary markers
	int numbdm = ptbdm.size();
	// number of segments
	int segmentadd = numseg;
	// a buffer vector
	vector<int> workVert(numseg,1);
	for (int i = 0; i < numseg; i = i + 1)
	{
		workVert[i] = i+1;
	}

	for (int k = 0; k < numseg; k = k + 1)
	{

		cout << "k="+k;

		int currentSeg = workVert[0];
		// first point in the current segment
		int point1 = seg2point[currentSeg-1][1];
		// second point in the current segment
		int point2 = seg2point[currentSeg-1][2];
		// segment boundary marker
		int bdm = seg2point[currentSeg-1][3];
		double x1; double x2; double y1; double y2;
		if (xcoord[point1-1] < xcoord[point2-1])
		{
			x1 = xcoord[point1 - 1];  y1 = ycoord[point1 - 1];
			x2 = xcoord[point2 - 1];  y2 = ycoord[point2 - 1];

		}
		else
		{
			 x1 = xcoord[point2 - 1];  y1 = ycoord[point2 - 1];
			 x2 = xcoord[point1 - 1];  y2 = ycoord[point1 - 1];
			 point1 = seg2point[currentSeg-1][2];
			 point2 = seg2point[currentSeg-1][1];
				
		}
		vector <int> positionx(2,0);
		
		// 2 for loops find out where the segment lies among the sectioniong planes defined in X
		for (int i = 0; i < Xlength; i = i + 1)
		{
			if (x1 <= X[i])
			{
				positionx[0] = i;
				break;
			}
			else if (x1 >= X[i] && i == (Xlength-1))
			{
				positionx[0] = i + 1;
			}
		}
		for (int i = 0; i < Xlength; i = i + 1)
		{
			if (x2 <= X[i])
			{
				positionx[1] = i;
				break;
			}
			else if (x2 >= X[i] && i == (Xlength-1))
			{
				positionx[1] = i + 1;
			}
		}

		sort(positionx.begin(), positionx.end());
		
		
		int deltax = abs(positionx[1] - positionx[0]);
		// if both points are in the same cell, then the segment is preserved in the x-direction
		if (deltax == 0)
		{
			workVert.push_back(currentSeg);
			workVert.erase(workVert.begin());
		}
		// otherwise, the segment must be split depending on the number of sections it passes through
		else
		{
			int newsegments = deltax + 1;
			int section = positionx[0];
			vector <double> coord1(2,0);
		    vector <double> coord2(2,0);
			for (int j = 0; j < newsegments; ++j)
			{
				//Added by AKANSHA
				if(j==3)
					cout << "break here";
				
				coord1.push_back(x1); coord1.push_back(y1);
				coord2.push_back(x2); coord2.push_back(y2);

				//Added by AKANSHA
				if(section >= Xlength-1)
					section=Xlength-1;

				double x3 = X[section];
				double y3 = Yintersect(x1,y1,x2,y2,x3);

				// replace original segment with shortened segment with same ID
				vector <double> newsegmentrow;
								
				//vector <double> newattrow(numatt-1,0);
				//Added by AKANSHA
				vector <double> newattrow(numatt,0);

			    if (j == 0)
				{
					// adding the new point with zero for attributes and boundary markers because it is
					// a created points
					xcoord.push_back(x3);
					ycoord.push_back(y3);
					ptbdm.push_back(0);
				    ptatt.push_back(newattrow);
					numpts = numpts + 1;

					newsegmentrow.push_back(currentSeg);
					newsegmentrow.push_back(point1);
					newsegmentrow.push_back(numpts);
					if (numbdm > 0)
					{
						newsegmentrow.push_back(bdm);
					}
					seg2point[currentSeg-1] = newsegmentrow;
					workVert.push_back(currentSeg);	
					segmentadd = segmentadd + 1;
				}

				// adds on middle segments
				else if (j > 0 && j < (newsegments - 1))
				{
					// adding the new point
					xcoord.push_back(x3);
					ycoord.push_back(y3);
					ptbdm.push_back(0);
					ptatt.push_back(newattrow);
					numpts = numpts + 1;

				    newsegmentrow.push_back(segmentadd);
					newsegmentrow.push_back(numpts-1);
					newsegmentrow.push_back(numpts);
					if (numbdm > 0)
					{
						newsegmentrow.push_back(bdm);
					}
					seg2point.push_back(newsegmentrow);
					workVert.push_back(segmentadd);
					segmentadd = segmentadd + 1;
				}
				// adds on the last segment
				else if ( j == (newsegments-1))
				{
					newsegmentrow.push_back(segmentadd);
					newsegmentrow.push_back(point2);
					newsegmentrow.push_back(numpts);
					newsegmentrow.push_back(bdm);
					
					seg2point.push_back(newsegmentrow);
					workVert.push_back(segmentadd);
					
					
				}
				section = section + 1;
				
				
			}
			workVert.erase(workVert.begin());
		}
	}
		
	numseg = segmentadd;
}

void SegmentSplitY(vector <double> & Y, vector< vector <double> > & seg2point, vector < double > & xcoord, vector <double> & ycoord, vector< vector <double> > & ptatt, vector<double> & ptbdm, int & numseg, int & numpts)
{
	// size of the Y vector
	int Ylength = Y.size();
	// number of attributes
	int numatt = ptatt[1].size();
	// number of boundary markers
	int numbdm = ptbdm.size();
	// number of segments
	int segmentadd = numseg;
	// a buffer vector
	vector<int> workHoriz(numseg,1);
	for (int i = 0; i < numseg; i = i + 1)
	{
		workHoriz[i] = i+1;
	}

	// takes each segment and sections it by horizontal section only
	for (int k = 0; k < numseg; ++k)
	{
		int currentSeg = workHoriz[0];
		// first point in the current segment
		int point1 = seg2point[currentSeg-1][1];
		// second point in the current segment
		int point2 = seg2point[currentSeg-1][2];
		// segment boundary marker
		int bdm = seg2point[currentSeg-1][3];
		double x1; double x2; double y1; double y2;
		if (ycoord[point1-1] < ycoord[point2-1])
		{
			x1 = xcoord[point1 - 1];  y1 = ycoord[point1 - 1];
			x2 = xcoord[point2 - 1];  y2 = ycoord[point2 - 1];
			point1 = seg2point[currentSeg-1][1];
			point2 = seg2point[currentSeg-1][2];
		}
		else
		{
			 x1 = xcoord[point2 - 1];  y1 = ycoord[point2 - 1];
			 x2 = xcoord[point1 - 1];  y2 = ycoord[point1 - 1];
			 point1 = seg2point[currentSeg-1][2];
			 point2 = seg2point[currentSeg-1][1];	
		}
		vector <int> positiony(2,0);

		// 2 for loops find out where the segment lies among the sectioning planes defined in Y
		for (int i = 0; i < Ylength; ++i)
		{
			if (y1 <= Y[i])
			{
				positiony[0] = i;
				break;
			}
			else if(y1 >= Y[i] && i == (Ylength - 1))
			{
				positiony[0] = i + 1;
			}
		}

		for (int i = 0; i < Ylength; ++i)
		{
			if (y2 <= Y[i])
			{
				positiony[1] = i;
				break;
			}
			else if (y2 >= Y[i] && i == (Ylength-1))
			{
				positiony[1] = i + 1;
			}
		}

		sort(positiony.begin(), positiony.end());
		
		int deltay = abs(positiony[1] - positiony[0]);
		// if both points are in the same cell, then the segment is preserved in the y-direction
		if (deltay == 0)
		{
			workHoriz.push_back(currentSeg);
			workHoriz.erase(workHoriz.begin());
		}
		// otherwise, the segment must be split depending on the number of sections it passes through
		else
		{
			int newsegments = deltay + 1;
			int section = positiony[0];
			vector <double> coord1(2,0);
			vector <double> coord2(2,0);
			for (int j = 0; j < newsegments; ++j)
			{
				coord1.push_back(x1); coord1.push_back(y1);
				coord2.push_back(x2); coord2.push_back(y2);

				//Added by AKANSHA
				if(section >= Ylength-1)
					section=Ylength-1;

				double y3 = Y[section];
				double x3 = Xintersect(x1,y1,x2,y2,y3);

				// replace original segment with shortened segment with same ID
				vector <double> newsegmentrow;
				
				//Updated by AKANSHA
				//vector <double> newattrow(numatt-1,0);
				vector <double> newattrow(numatt,0);

				if (j == 0)
				{
					// adding the new point with zero for attributes and boundary marker because it is
					// a created point
					xcoord.push_back(x3);
					ycoord.push_back(y3);
					ptbdm.push_back(0);
					ptatt.push_back(newattrow);
					numpts = numpts + 1;

					newsegmentrow.push_back(currentSeg);
					newsegmentrow.push_back(point1);
					newsegmentrow.push_back(numpts);
					if (numbdm > 0)
					{
						newsegmentrow.push_back(bdm);
					}
					seg2point[currentSeg-1] = newsegmentrow;
					workHoriz.push_back(currentSeg);
					segmentadd = segmentadd + 1;
				}

				// adds on middle segments
				else if (j > 0 && j < (newsegments-1))
				{
					// adding the new point
					xcoord.push_back(x3);
					ycoord.push_back(y3);
					ptbdm.push_back(0);
					ptatt.push_back(newattrow);
					numpts = numpts + 1;

					newsegmentrow.push_back(segmentadd);
					newsegmentrow.push_back(numpts-1);
					newsegmentrow.push_back(numpts);
					if (numbdm > 0)
					{
						newsegmentrow.push_back(bdm);
					}
					seg2point.push_back(newsegmentrow);
					workHoriz.push_back(segmentadd);
					segmentadd = segmentadd + 1;
				}

				// adds on last segment
				else if (j == (newsegments - 1))
				{
					
					newsegmentrow.push_back(segmentadd);
					newsegmentrow.push_back(point2);
					newsegmentrow.push_back(numpts);
					newsegmentrow.push_back(bdm);

					seg2point.push_back(newsegmentrow);
					workHoriz.push_back(segmentadd);
				}
				section = section + 1;
			}
			workHoriz.erase(workHoriz.begin());
		}
	}
	numseg = segmentadd;
}

// Produces the intersection points of the grid overlayed and specified by the user
void Grid(vector <double> & X, vector <double> & Y, vector <double> & xcoord, vector <double> & ycoord,
	vector < vector <double> > & ptatt, vector <double> & ptbdm, vector < vector <double> > & seg2point, 
	int & numatt, int & numpts, int & numbdm, int & numseg, int & numbdmseg)
{
	vector <double> Xoriginal = X;
	vector <double> Yoriginal = Y;

	double xmin = *std::min_element(xcoord.begin(), xcoord.end());
	double xmax = *std::max_element(xcoord.begin(), xcoord.end());
	double ymin = *std::min_element(ycoord.begin(), ycoord.end());
	double ymax = *std::max_element(ycoord.begin(), ycoord.end());

	// for the creation of the perimiter of the grid, we will need the minimum and maximum values for
	// each point
	X.push_back(xmin);
	X.push_back(xmax);
	Y.push_back(ymin);
	Y.push_back(ymax);

	sort(X.begin(), X.end());
	sort(Y.begin(), Y.end());

	int Xlength = X.size();
	int Ylength = Y.size();

	// stores new point IDs
	//Changed by Akansha
	//Previous: vector <int> ID ((Xlength*Ylength)-1, 0);
	vector <int> ID ((Xlength*Ylength), 0);

	// all grid point attributes will have a value of zero
	if (numatt > 0)
	{
		vector <double> newattrow(numatt-1,0);
	
		for (int i = 0; i < (Xlength*Ylength)-1; ++i)
		{
			ptatt.push_back(newattrow);
		}
	}
	// grid point boundary markers are -1
	if (numbdm > 0)
	{
		for (int i = 0; i < (Xlength*Ylength)-1; ++i)
		{
			ptbdm.push_back(-1);
		}
	}

	// adding in the coordinates of the gridpoints to xcoord and ycoord
	int gridcount = 1;
	int temp = numpts;
	for (int i = 0; i < Xlength; ++i)
	{
		for (int j = 0; j < Ylength; ++j)
		{
			xcoord.push_back(X[i]);
			ycoord.push_back(Y[j]);
			ID[gridcount-1] = numpts + gridcount;
			gridcount = gridcount + 1;
			temp = temp + 1;
		}
	}

	// now we add the segments of the grid
	gridcount = 0;
	int segcount = 1;
	int segID = numseg + 1;
	
	for (int i = 0; i < Xlength; ++i)
	{
		for (int j = 0; j < Ylength; ++j)
		{
			if ((j < Ylength - 1) || ( j == 0))
			{
				vector <double> newsegrow;
				newsegrow.push_back(segID);
				newsegrow.push_back(ID[gridcount]);
				newsegrow.push_back(ID[gridcount+1]);
				// we assign an grid boundary a value of -1.
				if (numbdmseg > 0)
				{
					newsegrow.push_back(-1);
				}
				seg2point.push_back(newsegrow);
				gridcount = gridcount + 1;
				segcount = segcount + 1;
				segID = segID + 1;
			}
			else
			{
				gridcount = gridcount + 1;
			}
		}
	}
	
	gridcount = 0;
	
	for (int i = 0; i < Ylength; ++i)
	{
		for (int j = 0; j < Xlength; ++j)
		{
			vector <double> newsegrow;
			if (j < (Xlength-1) || (j==0))
			{
				newsegrow.push_back(segID);
				newsegrow.push_back(ID[gridcount]);
				newsegrow.push_back(ID[gridcount + Ylength]);
				// we assign an grid boundary a value of -1.
				if (numbdmseg > 0)
				{
					newsegrow.push_back(-1);
				}
				seg2point.push_back(newsegrow);
				gridcount = gridcount + Ylength;
				segcount = segcount + 1;
				segID  = segID + 1;
			}
			else if (i < (Ylength-1))
			{
				gridcount = i+1;
			}
			else
			{
				break;
			}
		}
	}
	numpts = temp;
	numseg = segcount+numseg-1;
}
vector < vector < vector <double> > > PointSplit(vector <double> xcoord, vector <double> ycoord,
	vector <double> X, vector <double> Y, vector <vector <double> > ptatt, vector <double> ptbdm, 
	int numpts, int numatt, int numbdm)
{
	int Xlength = X.size();
	int Ylength = Y.size();
	vector < vector < vector <double> > > polycell; 
	for (int i = 0; i < (Xlength-1); ++ i)
	{
		for (int j = 0; j < (Ylength-1); ++j)
		{
			vector < vector <double> > minicell;
			for (int k = 0; k < numpts; ++k)
			{
				if ((xcoord[k]>=X[i]) && (xcoord[k]<=X[i+1]) && (ycoord[k]>=Y[j]) && (ycoord[k]<=Y[j+1]))
				{
					vector <double> cellrow;
					cellrow.push_back(k+1);
					cellrow.push_back(xcoord[k]);
					cellrow.push_back(ycoord[k]);
					if (numatt>0)
					{
						//vector <double> attrow  = ptatt[k];
						/*for (int w = 0; w < numatt; ++w)
						{
							cellrow.push_back(attrow[w]);
							cellrow.push_back(attrow[w]);
						}*/
					}
					if (numbdm>0)
					{
						cellrow.push_back(ptbdm[k]);
					}

					minicell.push_back(cellrow);
				}
				
			}
			polycell.push_back(minicell);
		}
	}
	return polycell;
}
vector < vector < vector <double> > > SegSplit(vector < vector < vector <double> > > polycell, 
	vector < vector <double> > seg2point, vector <double> X, vector <double> Y, int numbdm, int numseg)
{
	int Xlength = X.size();
	int Ylength = Y.size();
	vector < vector < vector <double> > > segcell;
	int count = 0;
	for (int i = 0; i < (Xlength-1); ++i)
	{
		vector <double> cellpoints;
		for (int j = 0; j < (Ylength-1); ++j)
		{
			vector < vector <double> > minicell;
			for (int a = 0; a < polycell[count].size(); ++a)
			{
				cellpoints.push_back(polycell[count][a][0]);
			}
			for (int k = 0; k < numseg; ++k)
			{
				vector <double> cellrow;
				int point1 = seg2point[k][1];
				int point2 = seg2point[k][2];
				int flag = 0;
				for (int a = 0; a < cellpoints.size(); ++a)
				{
					if (cellpoints[a] == point1)
					{
						flag = 1;
					}
					if (cellpoints[a] == point2)
					{
						flag = flag + 1;
					}
				}
				if (flag == 2)
				{
					cellrow.push_back(k+1);
					cellrow.push_back(point1);
					cellrow.push_back(point2);
					if (numbdm > 0)
					{
						cellrow.push_back(seg2point[k][3]);
					}
					minicell.push_back(cellrow);
				}
			}
			segcell.push_back(minicell);
			count = count + 1;
			cellpoints.clear();
		}	
		printVector(cellpoints);
	}
	return segcell;
}

// renumbers the points in each cell so that they correspond to the renumbered segment numbers
void Renumber (vector <double> & X, vector <double> & Y, vector < vector < vector <double> > > & polycell,
	vector < vector < vector <double> > > & segcell)
{
	int Xlength = X.size() - 2;
	int Ylength = Y.size() - 2;
	int numcells = (Xlength+1)*(Ylength+1);
	vector < vector <double> > minicellpoly;
	vector < vector <double> > minicellseg;
	for (int i = 0; i < numcells; ++i)
	{
		// the corresponding cells from polycell and segcell
		minicellpoly = polycell[i];
		minicellseg = segcell[i];
		// stores point ID's from each cell
		vector <double> points;
		// stores segments points from each cell
		vector < vector <double> > segmentpoints;
		for (int j = 0; j < minicellpoly.size(); ++j)
		{
			points.push_back(minicellpoly[j][0]);
		}
		// number of points in the cell
		int numpts = points.size();
		for (int j = 0; j < minicellseg.size(); ++j)
		{
			vector <double> segmentpointsrow;
			segmentpointsrow.push_back(minicellseg[j][1]);
			segmentpointsrow.push_back(minicellseg[j][2]);
			segmentpoints.push_back(segmentpointsrow);
		}

		for (int k = 0; k < numpts; ++k)
		{
			int coordinate = points[k];
			for (int a = 0; a < segmentpoints.size(); ++a)
			{
				for (int b = 0; b < segmentpoints[a].size(); ++b)
				{
					if (coordinate == segmentpoints[a][b])
					{
						segmentpoints[a][b] = k + 1;	
					}
				}
			}
			points[k] = k+1;
			minicellpoly[k][0] = points[k];
		}
		for (int a = 0; a < segmentpoints.size(); ++a)
		{
			minicellseg[a][1] = segmentpoints[a][0];
			minicellseg[a][2] = segmentpoints[a][1];
		}
		polycell[i] = minicellpoly;
		segcell[i] = minicellseg;
	}
}


	



int main(int argc, char** argv)
{

	std::cout<<"Entered main in polygridmain .";
	
	 string filename; 
	 vector <double> X;
     vector <double> Y;
	 
	 //cout << "Enter a filename: ";
	 //cin >> filename;
	 filename="simple_square.poly";

	 //cout << "Enter the number of x-partitions: ";
	 int Xlength;
	 //cin >> Xlength;
	 Xlength=3;

	 // creating the X partition vector
	 int counter;
	 counter = 0;
	 //while (counter<=(Xlength-1))
	 {
		//int xpartition;
        //Changed by AKANSHA
        double xpartition;
		//cout << "Enter a x-partition: ";
		//cin >> xpartition;
		//X.push_back(xpartition);
		X.push_back(0.2);
		X.push_back(0.5);
		X.push_back(0.9);
		counter = counter + 1;
	 }
	 sort(X.begin(), X.end());
	 

	 // creating the Y partition vector
	 //cout << "Enter the number of y-partitions: ";
	 int Ylength;
	 //cin >> Ylength;
	 Ylength=1;

	 counter = 0;
	 //while (counter<=(Ylength-1))
	 {
		//int ypartition;
        //Changed by AKANSHA
        double ypartition;
		//cout << "Enter a y-partition: ";
		//cin >> ypartition;
		//Y.push_back(ypartition);
		Y.push_back(0.5);
		counter = counter + 1;
	 }
	 
	 sort(Y.begin(),Y.end());

	  //open the input .poly file
	 ifstream inputfile(filename.c_str());

	 if(!inputfile.good())
	 {
		 stringstream str;
		 str<<"Unable to open file "<<filename<<".";
		 cout << str;
		 //throw runtime_error(str.str());
	 
	 }
	 //read data file
	 vector<double> v1;
	 while(!inputfile.eof())
	 {
		 string line;
		 getlinesc(inputfile,line);
		 trim(line);

		 string delimiter = " ";
		 vector <string> tokens;
		 Tokenize(line,tokens,delimiter);
		
		 if(!tokens.empty())
		 {
			 v1.reserve(tokens.size());
			 for(int i = 0; i < tokens.size(); ++i)
			 {
				 double element = String2Double(tokens[i]);
				 v1.push_back(element);
			 }
		 }
	 }
	 
	 inputfile.close();
	 
     // number of vertices in the input file
	 int numpts = v1[0];
	 // number of dimensions
	 int numdim = v1[1];
	 // number of attributes
	 int numatt = v1[2];
	 // number of boundary markers (either 0 or 1)
	 int numbdm = v1[3];
	 // number of elements in v1 that serve to define the file
	 int defline = 3;
	 // shift factor to get to each boundary marker stored in v1
	 int shiftbdm = 1 + numdim + numatt + numbdm;
	 // vector ptbdm stores all the boundary markers for each point
	 vector <double> ptbdm(numpts,0);
	 if (numbdm > 0)
	 {
		 for (int i = 0; i < numpts; i = i + 1)
		 {
			 ptbdm[i] = v1[defline + ((i+1)*shiftbdm)];
		 }
	 }
	 
	 
	 //vector ptatt stores all attributes for each point
     vector <double> att(numatt,0);
     vector< vector <double> > ptatt(numpts,att);
	 if (numatt > 0) 
	 {
		 // shift factor to get to the point attributes in v1
		 int shiftatt = 6;
		 for (int i = 0; i < numpts; i = i + 1)
		 {
			 
			 for (int j = 0; j < numatt; j = j + 1)
			 {
				 if (i == 0)
				 {
					ptatt[i][j] = v1[defline+4+j];
				 }
				 if (i > 0)
				 {
					 ptatt[i][j] = v1[defline + 4 + (i)*shiftatt + j-i];
					
				 }
			 }
		 }
	 }
	 	 

	 // vector xcoord stores all x-coordinates for all vertices in the input .poly file
	 vector <double> xcoord(numpts,0);
	 // vector ycoord stores all y-coordinates for all vertices in the input .poly file
	 vector <double> ycoord(numpts,0);
	 // shift factor to get to the next x or y coordinate in v1
	 int shiftpt = 3 + numatt + numbdm;

	 //Changed by Akansha
	 //Previous: for (int i = 0; i <= numpts; i = i + 1)
	 for (int i = 0; i < numpts; i = i + 1)
	 {
		 if (i == 0)
		 {
			 xcoord[i] = v1[defline+2];
			 ycoord[i] = v1[defline+3];
		 }
		 if(i > 0)
		 {
			 xcoord[i] = v1[defline+2+((i)*shiftpt)];
			 ycoord[i] = v1[defline+3+((i)*shiftpt)];
		 }
	 }

	 //getting to the second part of v1 that defines that segments
	 int shift2seg = defline+1+(numpts*shiftbdm); 
	 //  the second part of the .poly file that refers to the segments, holes, and regional attributes
	 vector <double> v2(v1.size()-shift2seg, 0);
	 int j = 0;
	 for (int i = 0; i < (v1.size()-shift2seg); i = i + 1)
	 {
		 v2[i] = v1[shift2seg+i];
		 j = j + 1;
	 }
	 

	 // number of segments in .poly file
	 int numseg = v2[0];
	 // number of boundary markers for each segments (0 or 1)
	 int numbdmseg = v2[1];
	 // number of elements in v2 that serve to define the segment portion of the .poly file
	 int deflineseg = 1;
	 // vector edbm defines the portion of the .poly file that refers to the segments
	 vector <double> edbm((deflineseg+((3+numbdmseg)*numseg)+1),0);
	 for (int i = 0; i <= (deflineseg+((3+numbdmseg)*numseg)); i = i + 1)
	 {
		 edbm[i] = v2[i];
	 }

	 
	 // shift factor to get to each point boundary marker stored in edbm
	 int shiftsegbdm = 3;
	 // vector of all boundary markers for all segments in the file
	 vector <double> bdmseg(numseg+1,0);
	 if (numbdmseg > 0)
	 {
		 for (int i = 0; i <= numseg; i = i + 1)
		 {
			 bdmseg[i] = edbm[deflineseg + i*(shiftsegbdm+1)];
		 }
	 }
	
	 // array that stores all information about segment point connections
	 vector <double> col(3+numbdmseg,0);
	 vector <vector <double> > seg2point(numseg,col);
	 // shift factor to get to vertices in edbm
	 int shiftsegpt = numbdmseg + 3;
	 for (int i = 0; i<numseg; i = i + 1)
	 {
		 for (int j = 0; j < 3+numbdmseg; j = j + 1)
		 {
			 if (i == 0)
			 {
				 seg2point[i][j] = edbm[deflineseg+1+j];
			 }
			 if (i>0)
			 {
				 seg2point[i][j] = edbm[(deflineseg+1) + ((i)*shiftsegpt) + j];
			 }
		 }
	 }

	 if (X.size() > 0)
	 {
		 // splitting the segments according to the user-defined X sections
		 SegmentSplitX(X,seg2point,xcoord,ycoord,ptatt,ptbdm,numseg,numpts);
	 }

	 if (Y.size()>0)
	 {
		 SegmentSplitY(Y,seg2point,xcoord,ycoord,ptatt,ptbdm,numseg,numpts);
	 }

     Grid(X,Y,xcoord,ycoord,ptatt,ptbdm,seg2point,numatt,numpts,numbdm,numseg,numbdmseg);

	 vector < vector < vector <double> > > polycell = PointSplit(xcoord,ycoord,X,Y,ptatt,ptbdm,numpts,numatt,numbdm);
	 vector < vector < vector <double> > > segcell = SegSplit(polycell,seg2point,X,Y,numbdm,numseg);

	 // renumbering the segments for each cell
	 int numcells = (Xlength+1)*(Ylength+1);
	 vector < vector <double> > minicell;
	 vector <double> cellrow;
	 for (int i = 0; i < numcells; ++i)
	 {
		 minicell = segcell[i];
		 for (int j = 0; j < minicell.size(); ++j)
		 {
			 minicell[j][0] = j+1;
		 }
		 segcell[i] = minicell;
		 minicell.clear();
	 }
	// portion that will renumber points for each cell corresponding to the segments
	Renumber(X,Y,polycell,segcell);
	
	// vector storing all filenames of the new filenames
	vector <string> filenames;
	for (int i = 0; i < numcells; ++i)
	{
		string tempfilename = filename;
		// converting i+1 to a string for the filename index
		string count;
		ostringstream convert;
		convert << i+1;
		count = convert.str();
		// adding count into the filename
		tempfilename.insert(filename.length()-5, count);
		// adding in the new filename into the filenames vector
		filenames.push_back(tempfilename);
	}

	FILE *outfile;
	
	// writing the new .poly files
	for (int i = 0; i < numcells; ++i)
	{
		ofstream newfile;
		string name = filenames[i];
		
		outfile = fopen(name.c_str(), "w+");
		  if (outfile == (FILE *) NULL) {
			printf("  Error:  Cannot create file %s.\n", name);
		  }
		
		
		//newfile.open(name);

		// the definition line for the poly section of the .poly file
		//ostringstream defline;
		//string deflinepoly; 
		//defline << polycell[i].size();
		//defline << 2;
		//defline << numatt;
		//defline << numbdm;
		//deflinepoly = defline.str();
		
		fprintf(outfile, "%d  %d  %d  %d\n", polycell[i].size(), 2,
          numatt, numbdm);

		//newfile << deflinepoly;
		//newfile << endl;

		for (int j = 0; j < polycell[i].size(); ++j)
		{
			string row;
			for (int k = 0; k < polycell[i][j].size(); ++k)
			{
				// taking in the elements corresponding to each file and converting them to strings in order
				// to be able to write them to the file
				string element;
				ostringstream convert;
				convert << polycell[i][j][k];
				element = convert.str();
				
				if (k < polycell[i][j].size() - 1)
				{
					row = row + element;
					row.push_back(' ');
				}
				else
				{
					row = row + element;
				}
			}
			fprintf(outfile, "%s\n", row);
			//newfile << row;
			//newfile<< endl;
		}
		//newfile.close();
		fclose(outfile);
	}

	 //for (int a = 0; a < polycell.size(); ++a)
	 //{
		// for (int b = 0; b < polycell[a].size(); ++b)
		// {
		//	 for (int c = 0; c < polycell[a][b].size(); ++c)
		//	 {
		//		 cout << polycell[a][b][c] << " ";
		//	 }
		//	 cout << endl;
		// }
		// cout << " " << endl;
	 //}

	 for (int a = 0; a < segcell.size(); ++a)
	 {
		 for (int b = 0; b < segcell[a].size(); ++b)
		 {
			 for (int c = 0; c < segcell[a][b].size(); ++c)
			 {
				 cout << segcell[a][b][c] << " ";
			 }
			 cout << endl;
		 }
		 cout << " " << endl;
	 }

	// cout << segcell.size() << endl;

	 

	return 0;
}


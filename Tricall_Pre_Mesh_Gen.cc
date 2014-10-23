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
#include "Tricall_Pre_Mesh_Gen.h"

////////////////////////////////////////////////////////////////////////////////
///
/// @brief Added by Akansha Kumar.
/// @This class performs the following functions, 1>Reads a PSLG file and 
/// @cutpoint details from the user, 2>Cuts the global segment into a set of
/// @specified local sub-segments,3>For each sub-segment data_transfer_object
/// @is populated, and 4>generateTriangleMesh function of tricall_utility is
/// @called.
///
/// @todo implement stitching in this class.
/// @todo implement multi-processors calls 
///
////////////////////////////////////////////////////////////////////////////////

Tricall_Pre_Mesh_Gen::Tricall_Pre_Mesh_Gen(void)
{
}


Tricall_Pre_Mesh_Gen::~Tricall_Pre_Mesh_Gen(void)
{
}


template <class T> void printVectorOfVectors(vector < vector < T > > v)
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

struct vector_accumulator
  : public stapl::p_object
{
private:
  double m_sum;
public:
  vector_accumulator(void)
    : m_sum(0.)
  {} 
  vector <vector <double> > recvBdrInfo;
  void receive_this(std::vector<double> const& v)
  { 
    recvBdrInfo.push_back(v);
  }
  vector <vector <double> > get_recvBdrInfo()
  {
  return recvBdrInfo;
  }

};

////////////////////////////////////////////////////////////////////////////////
///
/// @brief This method reads a PSLG file and cut point details
///  from the user and generates the PSLG files for each sub-segment.
///
///
/// @todo To test the code cutpoint details are hardcoded and are not accepted
///  by the user.
///
////////////////////////////////////////////////////////////////////////////////
int Tricall_Pre_Mesh_Gen::generateMesh()
{
  cout << "\nEnetered Tricall_Pre_Mesh_Gen generateMesh\n" ;
   
  int nprocs = stapl::get_num_locations();
  int myrank = stapl::get_location_id();
  
  
  Tricall_Utility tu;	
  Tricall_Pre_Mesh_Gen tripremeshgen;	
  
  string filename = "simple_square.poly"; 
  vector <double> X;
  X.push_back(.2); X.push_back(.5);
  sort(X.begin(), X.end());
  vector <double> Y;
  Y.push_back(.5); Y.push_back(.7);
  sort(Y.begin(), Y.end());
     
  //open the input .poly file
  ifstream inputfile(filename.c_str());
  
  if(!inputfile.good())
  {
   stringstream str;
   str<<"Unable to open file"<<filename<<".";
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

  for (int i = 0; i <= numpts; i = i + 1)
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
   // Creating the grid to ensure subdomain convexity
  Grid(X,Y,xcoord,ycoord,ptatt,ptbdm,seg2point,numatt,numbdm,numbdmseg);

  //X and Y include the minimum and maximum x and y values
  int subdomainrows = X.size()-1;
  int subdomaincols = Y.size()-1;

  // Vector that contains subdomain coordinates
  vector <string> subcoord;
  for (int i = 0; i < subdomainrows; ++i)
  {
    for (int j = 0; j < subdomaincols; ++j)
    {
      string combined;
      ostringstream merge;
      merge << i << ',' << j;
      combined = merge.str();
      subcoord.push_back(combined);
    }
  }

  // The global processor reference
  vector <int> GlobalProcessorReference;
  int proccount = 0;
  for (int i = 0; i < subcoord.size(); ++i)
  {
    if (proccount == nprocs)
    {
      proccount = 0;
    }
    GlobalProcessorReference.push_back(proccount);
    proccount = proccount + 1;
  }
  //A map linking subdomain integer ID to subdomain coordinate IDs
  map <string, int> SubCoordToInt;
  for (int i = 0; i < subcoord.size(); i=i+1)
  {
    string coord = subcoord[i];
    SubCoordToInt[coord] = i;
  }
  enum PossibleBoundaries {bottom=2,right=1,top=0,left=3};
    // A vector detailing which subdomains belong to a certain processor
  vector <int> rankvec;
  for (int i = 0; i < GlobalProcessorReference.size(); ++i)
  {
    if (GlobalProcessorReference[i] == myrank)
    {
      rankvec.push_back(i);
    }
  }
  

  vector_accumulator accumulator;
  vector <int> globalID;
  vector <vector <double> > pt2coord;
  map <int, vector<double> > SubDomainInfo;
  
  for (int i = 0; i < rankvec.size(); ++i)
  {	
    int CurrSubDom = rankvec[i];
    
    int numberofpoints;
    int numberofpointattributes;
    int numberofsegments;
    int numberofholes;
    int numberofregions;
    int pointCntr=0;
    char *filenamestr;
    int segCntr=0;
    int segMarkerCntr=0;
    
    //This is the maximum triangle area of a triangle element in 
    //the triangle mesh
    double maximumTriangleArea=0.02;
    
    //The switch code is temporary
    switch(i){
      case 0:maximumTriangleArea=0.002;
        break;
      case 1:maximumTriangleArea=0.004;
        break;
      case 2:maximumTriangleArea=0.008;
        break;
      case 3:maximumTriangleArea=0.016;
        break;
      case 4:maximumTriangleArea=0.032;
        break;
      case 5:maximumTriangleArea=0.016;
        break;
      case 6:maximumTriangleArea=0.004;
        break;
      case 7:maximumTriangleArea=0.002;
        break;
      case 8:maximumTriangleArea=0.004;
        break;
      case 9:maximumTriangleArea=0.008;
        break;
      case 10:maximumTriangleArea=0.016;
        break;
      case 11:maximumTriangleArea=0.032;
        break;
      case 12:maximumTriangleArea=0.016;
        break;
      case 13:maximumTriangleArea=0.004;
        break;
      case 14:maximumTriangleArea=0.002;
        break;	
    }
    
    vector <double> info;
    string sub = subcoord[CurrSubDom];
    int position = sub.find(',');
    // The x-coordinate where the subdomain is located on the grid
    string xloc = sub.substr(0,position);
    // The y-coordinate where the subdomain is located on the grid
    string yloc = sub.substr(position+1);
    //Converting these values to ints to use as pointers
    int x = atoi(xloc.c_str());
    int y = atoi(yloc.c_str());
    
    // The minimum X-value of the subdomain
    double Xmin = X[x];
    // The maximum X-value of the subdomain
    double Xmax = X[x+1];
    // The minimum Y-value of the subdomain
    double Ymin = Y[y];
    // The maximum Y-value of the subdomain
    double Ymax = Y[y+1];
    
    info.push_back(Xmin); info.push_back(Xmax);
    info.push_back(Ymin); info.push_back(Ymax);
    
    // We now explore the neighbors of the current subdomain
    int numneighbors;
    double rightbound=1;double leftbound=3;double topbound=0;double botbound=2;
    int RightNeighX=x+1; int RightNeighY=y;
    int LeftNeighX=x-1; int LeftNeighY =y;
    int TopNeighX=x; int TopNeighY=y+1;
    int BotNeighX=x; int BotNeighY=y-1;
    if (x == 0)
    {
      if (y == 0)
      {
        numneighbors = 2;
        int RightNeighID = Neigh(RightNeighX,RightNeighY,SubCoordToInt);
        int RightNeighProc = GlobalProcessorReference[RightNeighID];
        
        int TopNeighID = Neigh(TopNeighX,TopNeighY,SubCoordToInt);
        int TopNeighProc = GlobalProcessorReference[TopNeighID];
        
        info.push_back(numneighbors);
        info.push_back(TopNeighID); info.push_back(TopNeighProc);
        info.push_back(top);
        info.push_back(RightNeighID); info.push_back(RightNeighProc);
        info.push_back(right);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();			
      }
      else if (y>0 && y<(subdomaincols-1))
      {
        numneighbors = 3;
        int RightNeighID = Neigh(RightNeighX,RightNeighY,SubCoordToInt);
        int RightNeighProc = GlobalProcessorReference[RightNeighID];
        
        int TopNeighID = Neigh(TopNeighX,TopNeighY,SubCoordToInt);
        int TopNeighProc = GlobalProcessorReference[TopNeighID];
        
        int BotNeighID = Neigh(BotNeighX,BotNeighY,SubCoordToInt);
        int BotNeighProc = GlobalProcessorReference[BotNeighID];
        
        info.push_back(numneighbors);
        info.push_back(TopNeighID); info.push_back(TopNeighProc);
        info.push_back(top);
        info.push_back(RightNeighID); info.push_back(RightNeighProc);
        info.push_back(right);
        info.push_back(BotNeighID); info.push_back(BotNeighProc);
        info.push_back(bottom);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
      else
      {
        numneighbors = 2;
        int RightNeighID = Neigh(RightNeighX,RightNeighY,SubCoordToInt);
        int RightNeighProc = GlobalProcessorReference[RightNeighID];
        
        int BotNeighID = Neigh(BotNeighX,BotNeighY,SubCoordToInt);
        int BotNeighProc = GlobalProcessorReference[BotNeighID];

        info.push_back(numneighbors);
        info.push_back(RightNeighID); info.push_back(RightNeighProc);
        info.push_back(right);
        info.push_back(BotNeighID); info.push_back(BotNeighProc);
        info.push_back(bottom);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();				
      }
    }
    else if (x>0 && x<(subdomainrows-1))
    {
      if (y == 0)
      {
        numneighbors = 3;
        int RightNeighID = Neigh(RightNeighX,RightNeighY,SubCoordToInt);
        int RightNeighProc = GlobalProcessorReference[RightNeighID];
        
        int TopNeighID = Neigh(TopNeighX,TopNeighY,SubCoordToInt);
        int TopNeighProc = GlobalProcessorReference[TopNeighID];
        
        int LeftNeighID = Neigh(LeftNeighX,LeftNeighY,SubCoordToInt);
        int LeftNeighProc = GlobalProcessorReference[LeftNeighID];

        info.push_back(numneighbors);
        info.push_back(TopNeighID); info.push_back(TopNeighProc);
        info.push_back(top);
        info.push_back(RightNeighID); info.push_back(RightNeighProc);
        info.push_back(right);
        info.push_back(LeftNeighID); info.push_back(LeftNeighProc);
        info.push_back(left);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
      else if (y>0 && y<(subdomaincols-1))
      {
        numneighbors = 4;
        int RightNeighID = Neigh(RightNeighX,RightNeighY,SubCoordToInt);
        int RightNeighProc = GlobalProcessorReference[RightNeighID];
        
        int TopNeighID = Neigh(TopNeighX,TopNeighY,SubCoordToInt);
        int TopNeighProc = GlobalProcessorReference[TopNeighID];
        
        int LeftNeighID = Neigh(LeftNeighX,LeftNeighY,SubCoordToInt);
        int LeftNeighProc = GlobalProcessorReference[LeftNeighID];
        
        int BotNeighID = Neigh(BotNeighX,BotNeighY,SubCoordToInt);
        int BotNeighProc = GlobalProcessorReference[BotNeighID];
        
        info.push_back(numneighbors);
        info.push_back(TopNeighID); info.push_back(TopNeighProc);
        info.push_back(top);
        info.push_back(RightNeighID); info.push_back(RightNeighProc);
        info.push_back(right);
        info.push_back(BotNeighID); info.push_back(BotNeighProc);
        info.push_back(bottom);
        info.push_back(LeftNeighID); info.push_back(LeftNeighProc);
        info.push_back(left);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
      else
      {
        numneighbors = 3;
        int RightNeighID = Neigh(RightNeighX,RightNeighY,SubCoordToInt);
        int RightNeighProc = GlobalProcessorReference[RightNeighID];
        
        int LeftNeighID = Neigh(LeftNeighX,LeftNeighY,SubCoordToInt);
        int LeftNeighProc = GlobalProcessorReference[LeftNeighID];
        
        int BotNeighID = Neigh(BotNeighX,BotNeighY,SubCoordToInt);
        int BotNeighProc = GlobalProcessorReference[BotNeighID];
        
        info.push_back(numneighbors);
        info.push_back(RightNeighID); info.push_back(RightNeighProc);
        info.push_back(right);
        info.push_back(BotNeighID); info.push_back(BotNeighProc);
        info.push_back(bottom);
        info.push_back(LeftNeighID); info.push_back(LeftNeighProc);
        info.push_back(left);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
    }
    else
    {
      if (y == 0)
      {
        numneighbors = 2;
        int TopNeighID = Neigh(TopNeighX,TopNeighY,SubCoordToInt);
        int TopNeighProc = GlobalProcessorReference[TopNeighID];
        
        int LeftNeighID = Neigh(LeftNeighX,LeftNeighY,SubCoordToInt);
        int LeftNeighProc = GlobalProcessorReference[LeftNeighID];
        
        info.push_back(numneighbors);
        info.push_back(TopNeighID); info.push_back(TopNeighProc);
        info.push_back(top);
        info.push_back(LeftNeighID); info.push_back(LeftNeighProc);
        info.push_back(left);
                
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
      else if (y>0 && y<(subdomaincols-1))
      {
        numneighbors = 3;
        int TopNeighID = Neigh(TopNeighX,TopNeighY,SubCoordToInt);
        int TopNeighProc = GlobalProcessorReference[TopNeighID];
        
        int LeftNeighID = Neigh(LeftNeighX,LeftNeighY,SubCoordToInt);
        int LeftNeighProc = GlobalProcessorReference[LeftNeighID];
        
        int BotNeighID = Neigh(BotNeighX,BotNeighY,SubCoordToInt);
        int BotNeighProc = GlobalProcessorReference[BotNeighID];
        
        info.push_back(numneighbors);
        info.push_back(TopNeighID); info.push_back(TopNeighProc);
        info.push_back(top);
        info.push_back(BotNeighID); info.push_back(BotNeighProc);
        info.push_back(bottom);
        info.push_back(LeftNeighID); info.push_back(LeftNeighProc);
        info.push_back(left);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
      else
      {
        numneighbors = 2;
        int LeftNeighID = Neigh(LeftNeighX,LeftNeighY,SubCoordToInt);
        int LeftNeighProc = GlobalProcessorReference[LeftNeighID];
        
        int BotNeighID = Neigh(BotNeighX,BotNeighY,SubCoordToInt);
        int BotNeighProc = GlobalProcessorReference[BotNeighID];
        
        info.push_back(numneighbors);
        info.push_back(BotNeighID); info.push_back(BotNeighProc);
        info.push_back(bottom);
        info.push_back(LeftNeighID); info.push_back(LeftNeighProc);
        info.push_back(left);
        
        SubDomainInfo[CurrSubDom] =  info;
        info.clear();
      }
    }
    
    Data_Transfer_Object dto;
  
    stringstream conv;
    conv << CurrSubDom;
    string suffix = conv.str();
    string subfilename = filename;
    string fileprefix = filename.substr(0,filename.size()-5);		
    subfilename.insert(filename.size()-5,suffix);
    
    int length=strlen(subfilename.c_str());
    string canonical_name=subfilename.substr(0,length-5);		
    filenamestr=const_cast<char *>(canonical_name.c_str()); 
    // Necessary data containers for writing the .poly files
    vector <double> newxcoord;
    vector <double> newycoord;
    vector <vector<double> > newptatt;
    vector <double> newptbdm;
    vector <vector <double> > SubDomSegments;
    // The segments are split
    SegmentSplit(subfilename,Xmin,Xmax,Ymin,Ymax,seg2point,xcoord,ycoord,ptatt, ptbdm,numatt,numbdm,numbdmseg,newxcoord,newycoord,newptatt,newptbdm, SubDomSegments);
  
    vector <vector <double> > RegionalAttributes = RegionalAtt(fileprefix,Xmin,Xmax,Ymin,Ymax);
    
    // writing the new .poly files
    char *name = (char*)subfilename.c_str();
    ofstream newfile;
    
    newfile.open(name);
    // the definition for the node part of the .poly file
    stringstream deflinenode;
    string nodedefline;
    deflinenode << newxcoord.size() << " " << 2 << " " << numatt << " " << numbdm;
    nodedefline = deflinenode.str();
    newfile << nodedefline << endl;
    // The vertex information
    
    //populate local variables to be filled in the dto
    numberofpoints=newxcoord.size();
    numberofpointattributes=numatt;
    numberofsegments=0;
    numberofholes=0;
    numberofregions=0;
    numberofsegments=SubDomSegments.size();

    vector <double> trianglelist;
    vector <double> regionlist(numberofregions * 4,0);
    vector <double> pointlist;
    vector <double> pointattributelist(numberofpoints * numberofpointattributes,0);		
    vector <int> pointmarkerlist;					
    
    vector <int> segmentlist;
    vector <int> segmentmarkerlist;
    
    vector <double> edgelist;
    
    

    //Populate vertex details in the poly file
    //Populate vertex details in dto
    for (int j = 0; j < newxcoord.size(); ++j)
    {
      
      string row;
      stringstream conv;
      conv << j+1 << " " << newxcoord[j] << " " << newycoord[j];
      pointlist.push_back(newxcoord[j]);
      pointlist.push_back(newycoord[j]);
      
      /**
      if (numatt > 0)
      {
        for (int k = 0; k < numatt; ++j)
        {
          double att = newptatt[j][k];
          conv << " " << att;
        }
      }
      if (numbdm < 0)
      {
        double nodebdm = ptbdm[j];
        conv << " " << nodebdm;
      }
      **/
      
      //TODOD point marker
      pointmarkerlist.push_back(0);
      row = conv.str();
      newfile << row << endl;		
    }
    
    
    
    // The segment definition line
    stringstream defseg;
    string segdefline;
    defseg << SubDomSegments.size() << " " << numbdmseg;
    segdefline = defseg.str();
    newfile << segdefline << endl;
    for (int j = 0; j < SubDomSegments.size(); ++j)
    {
      string row;
      stringstream conv;
      for (int k = 0; k < SubDomSegments[j].size(); ++k)
      {
        conv << SubDomSegments[j][k];
        if (k < SubDomSegments[j].size()-1)
        {
          conv << " ";
          segmentlist.push_back(SubDomSegments[j][k]);
          segCntr++;
        }
        else
        {
           segmentmarkerlist.push_back(SubDomSegments[j][k]);
          segMarkerCntr++;
        }
      }
      row = conv.str();
      newfile << row << endl;
    }
    
    
    //Create and populate dto from local variables
    
    dto.numberofpoints=numberofpoints;
    dto.numberofpointattributes=numberofpointattributes;
    dto.pointlist=pointlist;
    dto.pointattributelist=pointattributelist;
    dto.pointmarkerlist=pointmarkerlist;
    dto.numberofholes=numberofholes;
    dto.numberofregions=numberofregions;
    dto.regionlist=regionlist;
    dto.maximumTriangleArea=maximumTriangleArea;
    
    dto.filenamestr=filenamestr;			
    //numberofsegments=0;
    dto.numberofsegments=numberofsegments;			
    dto.segmentlist=segmentlist;
    dto.segmentmarkerlist=segmentmarkerlist;	
    dto.logFlag=0;
    char *meshParam;
    string meshParamStr="enpraDPQ";
    meshParam=const_cast<char *>(meshParamStr.c_str()); 
    dto.meshCommandLine=meshParam;
    
    tu.generateTriangleMesh(dto);
    
    pointlist=dto.pointlist;
    edgelist=dto.edgelist;
    trianglelist=dto.trianglelist;
    int numberoftriangles=dto.numberoftriangles;
    int numberofcorners=dto.numberofcorners;
    int updatednumberoftriangles=numberoftriangles;
    numberofcorners=3;
    
    int numberofnodes=pointlist.size();
    int updatednumberofnodes=numberofnodes;
    
    
    //triangle_ptr=array of triangle elements
    //Each element contains numberofcorners node indices and a flag(valid - 1 or invalid - 0)
    int *triangle_ptr;
    triangle_ptr = (int *) malloc(sizeof(int) * (numberoftriangles*(numberofcorners+1)));	
    for(int iTrianIter=0; iTrianIter < numberoftriangles;iTrianIter++)
    {
      (*(triangle_ptr + (iTrianIter*(numberofcorners+1)) + 0)) = trianglelist[(iTrianIter*(numberofcorners)) + 0];
      (*(triangle_ptr + (iTrianIter*(numberofcorners+1)) + 1)) = trianglelist[(iTrianIter*(numberofcorners)) + 1];
      (*(triangle_ptr + (iTrianIter*(numberofcorners+1)) + 2)) = trianglelist[(iTrianIter*(numberofcorners)) + 2];
      (*(triangle_ptr + (iTrianIter*(numberofcorners+1)) + 3)) = 1;
    }
    
    
    //This is the index added to the existing nodes to obtain the index of a new node
    int newNodeIndexAdder=0;
    
    vector_accumulator accumulator;
    
    //Now read the neighbouring details of the segment
    map<int, vector <double> >::iterator it;
    it = SubDomainInfo.find(CurrSubDom);
    vector <double> neigh_segment_proc_vec1=it->second;
    int number_of_neighbors=neigh_segment_proc_vec1[4];
    for(int iNeighIter=0; iNeighIter < number_of_neighbors; iNeighIter++ )
    {
      int neigh_seg=neigh_segment_proc_vec1[5 + (iNeighIter*3)];
      int neigh_seg_proc=neigh_segment_proc_vec1[6 + (iNeighIter*3)];
      int neigh_location=neigh_segment_proc_vec1[7 + (iNeighIter*3)];
              
      //Get the boundary elements of current segment and processor
      //location_of_boundary:: (0,top),(1,right),(2,bottom),(3,left)
      double bdrVal=0;
      if(neigh_location == top)
        bdrVal=neigh_segment_proc_vec1[3];
      if(neigh_location == right)
        bdrVal=neigh_segment_proc_vec1[1];
      if(neigh_location == bottom)
        bdrVal=neigh_segment_proc_vec1[2];
      if(neigh_location == left)
        bdrVal=neigh_segment_proc_vec1[0];
      vector<double> boundPtList;	
      boundPtList.push_back(myrank);
      boundPtList.push_back(CurrSubDom);
      boundPtList.push_back(neigh_location);
      
      //This holds the  
      vector <int> orgBoundNodeIndex;
      
      //getBoundaryNodes() returns, the x_values on the boundary if the boundary on the top or bottom, 
      //the y-values on the boundary if the boundary is on the left or right, and
      //other attributes
      Tricall_Pre_Mesh_Gen tricall;
      tricall.getBoundaryNodes(pointlist,neigh_location,bdrVal,boundPtList,myrank,CurrSubDom,neigh_seg_proc,neigh_seg,numberofnodes,orgBoundNodeIndex);
      
      //Send values to another processor having a neighboring segment
      //TAG:- Current Segment Id
      //DEST:- Destination processor
      stapl::async_rmi(neigh_seg_proc,accumulator.get_rmi_handle(), &vector_accumulator::receive_this,boundPtList);
            
    }
    
    stapl::rmi_fence();
    
    //Recieving
    vector< vector<double> > rcvdBdrInfo = accumulator.get_recvBdrInfo();
    //printVectorOfVectors(rcvdBdrInfo);
    
    //Receiving data from neighboring processors
    for(int iNeighIter=0; iNeighIter < rcvdBdrInfo.size(); iNeighIter++ )
    {
      int neigh_seg=rcvdBdrInfo.at(iNeighIter).at(1);
      int neigh_seg_proc=rcvdBdrInfo.at(iNeighIter).at(0);
      int neigh_location=-2;
      
      for(int iNeighIter=0; iNeighIter < number_of_neighbors; iNeighIter++ )
      {
        if(neigh_seg == neigh_segment_proc_vec1[5 + (iNeighIter*3)] && neigh_seg_proc == neigh_segment_proc_vec1[6 + (iNeighIter*3)])
        {
          neigh_location=neigh_segment_proc_vec1[7 + (iNeighIter*3)];
        }
      }
              
      //Declare a vector to store data obtained after MPI Recv
      vector<double> neigh_point_list=rcvdBdrInfo.at(iNeighIter);
      int number_amount=neigh_point_list.size()-3;
      
      
      //Get the boundary elements of current segment and processor
      //location_of_boundary:: (0,top),(1,right),(2,bottom),(3,left)
      double bdrVal=0;
      if(neigh_location == top)
        bdrVal=neigh_segment_proc_vec1[3];
      if(neigh_location == right)
        bdrVal=neigh_segment_proc_vec1[1];
      if(neigh_location == bottom)
        bdrVal=neigh_segment_proc_vec1[2];
      if(neigh_location == left)
        bdrVal=neigh_segment_proc_vec1[0];
      vector<double> boundPtList;	
      
      //This holds the  
      vector <int> orgBoundNodeIndex;
      
      //Key is the co-ordinate and value is the index
      map<int,double> nodesDataMap;
      
      //getBoundaryNodes() returns, the x_values on the boundary if the boundary on the top or bottom, 
      //the y-values on the boundary if the boundary is on the left or right, and
      //other attributes
      getBoundaryNodes(pointlist,neigh_location,bdrVal,boundPtList,myrank,CurrSubDom,neigh_seg_proc,neigh_seg,numberofnodes,orgBoundNodeIndex);
      
      updateContainers(nodesDataMap,boundPtList,orgBoundNodeIndex);
      
      vector<int>::iterator orgBoundNodeIter1;
      
      //Now parse through elements and obtain the boundary elements
      //Size of bndrElementsIndices is the size of boundary elements 
      vector<int> bndrElementsIndices;				
      //Value=element orig index,index of node1 on boundary,index of node2 on boundary, and index of the third node
      vector<int> bdrElementsVec;
      
      for(int iElementIter=0;iElementIter < numberoftriangles;iElementIter++)
      {
      
        if((*(triangle_ptr + (iElementIter*(numberofcorners+1)) + 3)) == 1)
        {
          //The triangle is a valid triangle.
          vector<int> indicesVec;
          int thirdNode=-1;
        
          int node1=(*(triangle_ptr + (iElementIter*(numberofcorners+1)) + 0));
          int node2=(*(triangle_ptr + (iElementIter*(numberofcorners+1)) + 1));
          int node3=(*(triangle_ptr + (iElementIter*(numberofcorners+1)) + 2));
          
          int cntr=0;
          
          orgBoundNodeIter1=find(orgBoundNodeIndex.begin(),orgBoundNodeIndex.end(),node1);
          if(orgBoundNodeIter1!=orgBoundNodeIndex.end())
          {
            cntr++;
            indicesVec.push_back(node1);
          }
          else
          {
            thirdNode=node1;
          }
          
          orgBoundNodeIter1=find(orgBoundNodeIndex.begin(),orgBoundNodeIndex.end(),node2);
          if(orgBoundNodeIter1!=orgBoundNodeIndex.end())
          {
            cntr++;	
            indicesVec.push_back(node2);
          }
          else
          {
            thirdNode=node2;
          }
          
          orgBoundNodeIter1=find(orgBoundNodeIndex.begin(),orgBoundNodeIndex.end(),node3);
          if(orgBoundNodeIter1!=orgBoundNodeIndex.end())
          {
            cntr++;
            indicesVec.push_back(node3);
          }
          else
          {
            thirdNode=node3;
          }

          if(cntr == 2)
          {
            indicesVec.push_back(thirdNode);
            bndrElementsIndices.push_back(iElementIter);
            
            bdrElementsVec.push_back(iElementIter);
            bdrElementsVec.push_back(indicesVec[0]);
            bdrElementsVec.push_back(indicesVec[1]);
            bdrElementsVec.push_back(indicesVec[2]);
            
          }
        }	
      }
      
      
      //Start stitching				
      double thresHold=1.0E-05;
      //Parse through all boundary elements of neighboring elements
      map<int,double >::iterator nodeIt;
      //neigh_point_list.size()
      for(int iNeighIter=0;iNeighIter < neigh_point_list.size();iNeighIter++)
      {
        double value=neigh_point_list[iNeighIter];
        
        //Parse through each boundary element
        int iBdrElemIter=0;
        while(iBdrElemIter < (bdrElementsVec.size()/4))			
        {
          //Key is the element index	
          int key = bdrElementsVec[(iBdrElemIter*4)+0];
          //Value contain indices of the two boundary nodes
          vector<int> map_value;
          map_value.push_back(bdrElementsVec[(iBdrElemIter*4)+1]);
          map_value.push_back(bdrElementsVec[(iBdrElemIter*4)+2]);
          map_value.push_back(bdrElementsVec[(iBdrElemIter*4)+3]);
          
          nodeIt=nodesDataMap.find(map_value[0]);
          double nodeVal1=nodeIt->second;
          nodeIt=nodesDataMap.find(map_value[1]);
          double nodeVal2=nodeIt->second;
          double tempVal=0;
          
          
          
          
          //cout << key << "###" << value << "," << map_value[0] << "," << map_value[1] << "::" << nodeVal1 << "," << nodeVal2 << "\n";	
          
          
          
          if( ( (nodeVal1 > nodeVal2) && ((value > nodeVal2) && (value < nodeVal1) && ((value - nodeVal2) > thresHold) && ((nodeVal1 - value) > thresHold))   )  ||  
          ( ((nodeVal1 <= nodeVal2)) && ((value > nodeVal1) && (value < nodeVal2) && ((value - nodeVal1) > thresHold) && ((nodeVal2 - value) > thresHold))   ) )
          {
          
            double xCoord=0.0;
            double yCoord=0.0;
            if(neigh_location ==top || neigh_location==bottom)
            {
              xCoord=value;
              yCoord=bdrVal;
            }
            if(neigh_location ==right || neigh_location==left)
            {
              xCoord=bdrVal;
              yCoord=value;
            }
            
            
            if(neigh_location ==top || neigh_location==bottom||neigh_location ==right || neigh_location==left)
            {
              //Add a node
              cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << "Adding a node" << "\n";
              
              //Resize and update point_ptr 
              //Resizing is done using a temporary pointer array
              pointlist.push_back(xCoord);
              pointlist.push_back(yCoord);
              updatednumberofnodes=updatednumberofnodes+2;
              
              newNodeIndexAdder++;
            }
          
          }
          iBdrElemIter++;
        }
        
      }
    }
    
    //Start
    //Now call the triangle library again with the same boundary nodes
    //and the same maximum area
    //Create updated triangles - Start
    Tricall_Utility tu1;
    Data_Transfer_Object dto1;
    cout <<"-----------------------------------------------------------------------------------------------------------------------" << "\n" ;
    int numberofpoints1=updatednumberofnodes/2;
    int numberofpointattributes1=0;
    int numberofsegments1;
    int numberofholes1=0;
    int numberofregions1=0;
    
    vector <double> trianglelist1;
    vector <double> regionlist1(numberofregions1 * 4,0);
    vector <double> pointlist111(numberofpoints1 * 2,0);
    vector <double> pointattributelist1(numberofpoints1 * numberofpointattributes1,0);		
    vector <int> pointmarkerlist1(numberofpoints1 ,0);		
    vector <int> segmentlist1(numberofsegments1 * 2 ,0);
    vector <int> segmentmarkerlist1(numberofsegments1 ,0);
    vector <double> edgelist1;
    //Populate vertex details in the poly file
    //Populate vertex details in dto
    for (int j = 0; j < (numberofpoints1*2); j++)
    {
      pointlist111[j]=pointlist[j];
    }						
    //Create and populate dto from local variables			
    dto1.numberofpoints=numberofpoints1;
    dto1.numberofpointattributes=numberofpointattributes1;
    dto1.pointlist=pointlist111;
    dto1.pointattributelist=pointattributelist1;
    dto1.pointmarkerlist=pointmarkerlist1;
    dto1.numberofholes=numberofholes1;
    dto1.numberofregions=numberofregions1;
    dto1.regionlist=regionlist1;
    dto1.maximumTriangleArea=maximumTriangleArea;
    canonical_name.append("_st1");
    char *filenamestr1;
    filenamestr1=const_cast<char *>(canonical_name.c_str()); 
    dto1.filenamestr=filenamestr1;		
    dto1.numberofsegments=numberofsegments1;		
    dto1.segmentlist=segmentlist1;	
    dto1.segmentmarkerlist=segmentmarkerlist1;
    dto1.edgelist=edgelist1;
    dto1.logFlag=999;
    char *meshParam1;			
    string meshParamStr1="a";
    stringstream numStrStream;
    numStrStream << maximumTriangleArea;
    meshParamStr1.append(numStrStream.str());
    meshParamStr1.append("eAnDYPQ");
    meshParam1=const_cast<char *>(meshParamStr1.c_str()); 
    dto1.meshCommandLine=meshParam1;
    
    tu1.generateMeshWithFixedBoundaries(dto1);
    //End
      
      
    
    
    
     /**for (int a = 0; a < segcell.size(); ++a)
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
     }**/

    
    
  }
  
  
}

void Tricall_Pre_Mesh_Gen::updateContainers(map<int,double>& nodesDataMap,vector <double>& boundPtList,vector <int>& orgBoundNodeIndex)
{
  for (int i = 0; i < (boundPtList.size()); i++)
  {
    nodesDataMap[orgBoundNodeIndex[i]]=boundPtList[i];
    
  }	
  
}

 //bdrIndex=0(top),1(right),2(bottom),3(left) -input
 //bdrVal = value on the boundary - input
 //number_count = total number of points.
 //boundPtList containes bundary points (x-coordinates - (top and bottom), y-coordinates - (left and right)) - output
 //orgBoundNodeIndex: containes indices of original boundary nodes - output
 //points_data_frm_proc contains all the nodes as input
 //TODO::int myrank,int i,int neigh_seg_proc,int neigh_seg :: Just for testing .. will be removed
 void Tricall_Pre_Mesh_Gen::getBoundaryNodes(vector<double> &points_data_frm_proc , int bdrIndex, double bdrVal, vector<double> &boundPtList, int my_id,int i,int neigh_seg_proc,int neigh_seg,int number_count,vector<int> &orgBoundNodeIndex)
 {
  double epsilon=1.0E-5;
  //cout << my_id << i << neigh_seg_proc << neigh_seg << "Started getBoundaryNodes method" << "\n";
  
  if(bdrIndex == 0 || bdrIndex == 2)
  {
    for(int iPoint=0; iPoint < (number_count/2) ; iPoint++)
    {
      double y_plus_bdr=bdrVal;
      double x_val=points_data_frm_proc[iPoint*2];
      double y_val=points_data_frm_proc[(iPoint*2) + 1];
      //cout << x_val << "~~" << y_val << "\n";
      if( y_val >= (y_plus_bdr - epsilon) && y_val <= (y_plus_bdr + epsilon))
      {
        boundPtList.push_back(x_val);
        orgBoundNodeIndex.push_back(iPoint);
        //boundPtList.push_back(y_val);
        //cout << x_val << "::" << y_val << "\n";
      }
    }
  }
  if(bdrIndex == 1 || bdrIndex == 3)
  {
    for(int iPoint=0; iPoint < (number_count/2) ; iPoint++)
    {
      double x_plus_bdr=bdrVal;
      double x_val=points_data_frm_proc[iPoint*2];
      double y_val=points_data_frm_proc[(iPoint*2) + 1];
      //cout << x_val << "~~" << y_val << "\n";
      if( x_val >= (x_plus_bdr - epsilon) && x_val <= (x_plus_bdr + epsilon))
      {
        //boundPtList.push_back(x_val);
        boundPtList.push_back(y_val);
        orgBoundNodeIndex.push_back(iPoint);
        //cout << x_val << "::" << y_val << "\n";
      }
    }
  }
  //cout << my_id << i << neigh_seg_proc << neigh_seg << "Ended getBoundaryNodes method" << "\n";
 }
 
 

// Utility functions used for testing purposes

inline void Tricall_Pre_Mesh_Gen::getlinesc(std::ifstream& inputfile, std::string& line)
{
  std::getline(inputfile, line);
  if(line[0] == COMMENT_CHAR)
  getlinesc(inputfile, line);
}
inline double Tricall_Pre_Mesh_Gen::String2Double(std::string s)
{
  double val;
  std::istringstream iss(s);
  iss >> val;
  return val;
}
inline void Tricall_Pre_Mesh_Gen::Tokenize(std::string& line, std::vector<std::string>& tokens,
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
void Tricall_Pre_Mesh_Gen::trim(string&s)
{
  size_t p = s.find_first_not_of("\t");
  s.erase(0,p);

  p = s.find_last_not_of("\t");
  if(string::npos != p)
    s.erase(p+1);
}
template <class T> void printVector(vector <T> v)
{
  for (int i = 0; i < v.size(); i = i + 1)
  {
    cout << "Element " << i << " = " << v[i] << endl;
  }
}


// Here are functions used within the operator code

void Tricall_Pre_Mesh_Gen::load_geom(int & ntri, vector <double> & AttID, vector < vector <int> > & tri2edge,vector < vector <int> > & edge2tri, vector < vector <int> > & edge2pt, vector < vector <int> > & tri2pt,vector < vector <double> > & pt2coord, string & fileprefix, int & index)
{
  string newindex;
  ostringstream convert;
  convert << index;
  newindex = convert.str();
  string filenameedge = fileprefix+"."+newindex+".edge";
  string filenamenode = fileprefix+"."+newindex+".node";
  string filenameele = fileprefix+"."+newindex+".ele";

  // opening the edge file and reading in the values
  ifstream edgefile(filenameedge.c_str());
  if (!edgefile.good())
  {
    stringstream str;
    str << "Unable to open file"<<filenameedge<<".";
  }
  // read in the edge file
  vector <double> v1;
  while(!edgefile.eof())
  {
    string line;
    getlinesc(edgefile,line);
    trim(line);
    string delimiter = " ";
    vector <string> tokens;
    Tokenize(line,tokens,delimiter);
    if(!tokens.empty())
    {
      v1.reserve(tokens.size());
      for (int i = 0; i < tokens.size(); ++i)
      {
        double element = String2Double(tokens[i]);
        v1.push_back(element);
      }
    }
  }
  edgefile.close();
  // the number of edges
  int nedge = v1[0];
  // the boundary markers
  int bm = v1[1];
  // edge numer array
  vector <int> ednbr;
  for(int i = 0; i < nedge; ++i)
  {
    ednbr.push_back(i+1);
  }
  // the edge to point array
  vector < vector <int> > e2p; 
  // edge boundary marker array
  vector <double> edbmar;
  int k = 2;
  int del = 3+bm;
  for (int j = 0; j < nedge; ++j)
  {
    int shift = k+(j*del);
    int point1 = v1[shift+1];
    int point2 = v1[shift+2];
    if (bm>0)
    {
      int bdm = v1[shift+3];
      edbmar.push_back(bdm);
    }	
    vector <int> newrow;
    newrow.push_back(point1);
    newrow.push_back(point2);
    e2p.push_back(newrow);
  }

  // opening the node file and read in the value
  ifstream nodefile(filenamenode.c_str());
  if (!nodefile.good())
  {
    stringstream str;
    str << "Unable to open file"<<filenamenode<<".";
  }
  // read in the node file
  vector <double> v2;
  while(!nodefile.eof())
  {
    string line;
    getlinesc(nodefile,line);
    trim(line);
    string delimiter = " ";
    vector <string> tokens;
    Tokenize(line,tokens,delimiter);
    if(!tokens.empty())
    {
      v2.reserve(tokens.size());
      for (int i = 0; i < tokens.size(); ++i)
      {
        double element = String2Double(tokens[i]);
        v2.push_back(element);
      }
    }
  }
  nodefile.close();

  // the number of points
  int npts = v2[0];
  // the number of dimensions
  int dim = v2[1];
  // the number of attributes
  int att = v2[2];
  // the number of boundary markers
  bm = v2[3];
  // the point number array
  vector <int> ptnbr;
  for (int i = 0; i < npts; ++i)
  {
    ptnbr.push_back(i+1);
  }
  // the coordinates of each point
  vector <vector <double> > p2c;
  // the attributes of each point
  vector <vector <double> > patt;
  // the boundary marker of each point
  vector <double> pbdmar;

  k = 4;
  del = 3+att+bm;
  for (int j = 0; j < npts; ++j)
  {
    int shift = k+(j*del);
    double xcoord = v2[shift+1];
    double ycoord = v2[shift+2];
    vector <double> newrow;
    newrow.push_back(xcoord);
    newrow.push_back(ycoord);
    p2c.push_back(newrow);
    vector <double> attrow;
    if (att > 0)
    {
      for (int k = 0; k < att; ++k)
      {
        attrow.push_back(v2[shift+3+k]);
      }
      patt.push_back(attrow);
    }
    if (bm>0)
    {
      pbdmar.push_back(v2[shift+2+att+1]);
    }
  }

  // open the element file and read in the values
  ifstream elefile(filenameele.c_str());
  if (!elefile.good())
  {
    stringstream str;
    str << "Unable to open file"<<filenameele<<".";
  }
  // read in the element file
  vector <double> v3;
  while(!elefile.eof())
  {
    string line;
    getlinesc(elefile,line);
    trim(line);
    string delimiter = " ";
    vector <string> tokens;
    Tokenize(line,tokens,delimiter);
    if(!tokens.empty())
    {
      v3.reserve(tokens.size());
      for (int i = 0; i < tokens.size(); ++i)
      {
        double element = String2Double(tokens[i]);
        v3.push_back(element);
      }
    }
  }
  elefile.close();

  // the number of triangles
  ntri = v3[0];
  // the number of dimensions
  dim = v3[1];
  // the number of attributes
  att = v3[2];
  // the triangle number array
  vector <int> trinbr;
  for (int i = 0; i < ntri; ++i)
  {
    trinbr.push_back(i+1);
  }
  // the triangle to point array
  vector < vector <int> > t2p;
  // the triangle attribute array
  vector <double> tatt;

  k = 3;
  del = 4+att;
  for (int j = 0; j < ntri; ++j)
  {
    int shift = k+(j*del);
    int tri1 = v3[shift+1];
    int tri2 = v3[shift+2];
    int tri3 = v3[shift+3];
    if (att>0)
    {
      tatt.push_back(v3[shift+3+att]);
    }
    vector <int> newrow;
    newrow.push_back(tri1);
    newrow.push_back(tri2);
    newrow.push_back(tri3);
    t2p.push_back(newrow);
  }
  AttID = tatt;
  edge2pt = e2p;
  tri2pt = t2p;
  pt2coord = p2c;
  vector <int> newedgerow;
  vector <int> row1;
  vector <int> row2;
  vector <int> indrow;
  for (int i = 0; i < nedge; ++i)
  {	
    int point1 = edge2pt[i][0];
    int point2 = edge2pt[i][1];
    // finding on what rows of t2p the two points exist on
    for (int j = 0; j < tri2pt.size(); ++j)
    {
      for (int k = 0; k < tri2pt[j].size(); ++k)
      {
        if (tri2pt[j][k] == point1)
        {
          row1.push_back(j);
        }
        if (tri2pt[j][k] == point2)
        {
          row2.push_back(j);
        }
      }
    }
    for (int a = 0; a < row1.size(); ++a)
    {
      for (int b = 0; b < row2.size(); ++b)
      {
        if (row1[a] == row2[b])
        {	
          int tri = row1[a];
          indrow.push_back(tri);
        }
      }
    }
    for (int k = 0; k < indrow.size(); ++k)
    {
      int tri = indrow[k]+1;
      newedgerow.push_back(tri);		
    }
    if (indrow.size() == 1)
    {
      newedgerow.push_back(0);
    }
    edge2tri.push_back(newedgerow);
    newedgerow.clear();
    row1.clear(); row2.clear();
    indrow.clear(); 
  }

  //creating tri2egde

  for (int i = 0; i < ntri; ++i)
  {
    for (int j = 0; j < edge2tri.size(); ++j)
    {
      int tri1 = edge2tri[j][0];
      int tri2 = edge2tri[j][1];
      if (tri1 == i+1)
      {
        row1.push_back(j+1);
      }
      if (tri2 == i+1)
      {
        row1.push_back(j+1);
      }
    }
    tri2edge.push_back(row1);
    row1.clear();
  }
}

// Given two points that form a line (input as a 1x2 [x y] vector) and a y-plane that intersects the segment formed by the two points, Xintersect finds the x-coordinate where the segment and plane intersect
double Tricall_Pre_Mesh_Gen::Xintersect(double x1, double y1, double x2, double y2, double y3)
{
  
  // slope of the segment
  double m = (y2-y1)/(x2-x1);

  double x3 = ((y3-y1)/m) + x1;

  return x3;
}
// Given that two points form a line, and a x-plane (x3) intersects that segment,
// this function finds the y coordinate where the line and plane intersect
double Tricall_Pre_Mesh_Gen::Yintersect(double x1, double y1, double x2, double y2, double x3)
{
  
  //slope of the segment
  double m = (y2-y1)/(x2-x1);

  double y3 = (m*(x3-x1)) + y1;

  return y3;
}

vector <vector<double> > Tricall_Pre_Mesh_Gen::SortSegments(vector <vector<double> > potentials, vector <double> ycoord)
{
  vector <double> newline;
  vector <vector<double> > pointmarkers;
  for (int i = 0; i < potentials.size(); ++i)
  {
    int point1 = potentials[i][1];
    int point2 = potentials[i][2];
    double y1 = ycoord[point1-1];
    double y2 = ycoord[point2-1];
    int m1; int m2;
    if (y1 < y2)
    {
      m1 = 1;
      m2 = 2;
    }
    else
    {
      m1 = 2;
      m2 = 1;
    }
    newline.push_back(y1);
    newline.push_back(point1);
    newline.push_back(m1);
    pointmarkers.push_back(newline);
    newline.clear();
    newline.push_back(y2);
    newline.push_back(point2);
    newline.push_back(m2);
    pointmarkers.push_back(newline);
    newline.clear();
  }
  sort(pointmarkers.begin(),pointmarkers.end());
  return pointmarkers;
}

void Tricall_Pre_Mesh_Gen::Grid(vector <double> & X, vector <double> & Y, vector <double> & xcoord, vector <double> & ycoord,vector <vector<double> > & ptatt, vector <double> & ptbdm, vector <vector <double> > & seg2point, int & numatt, int & numbdm,int & numbdmseg)
{
  
  double xmin = *std::min_element(xcoord.begin(), xcoord.end());
  double xmax = *std::max_element(xcoord.begin(), xcoord.end());
  double ymin = *std::min_element(ycoord.begin(), ycoord.end());
  double ymax = *std::max_element(ycoord.begin(), ycoord.end());

  // for the creation of the perimeter of the grid, we will need the minimum and maximum values for
  // each point
  X.push_back(xmin);
  X.push_back(xmax);
  Y.push_back(ymin);
  Y.push_back(ymax);
  
  sort(X.begin(), X.end());
  sort(Y.begin(), Y.end());
  
  // We know check if any portion of grid segments along the X-cut planes exist, because we must preserve their boundary values.
  vector <double> newseg;
  vector <vector<double> > potentials;
  vector <double> newattrow;
  for (int i = 0; i < X.size(); ++i)
  {
    for (int j = 0; j < seg2point.size(); ++j)
    {
      int point1 = seg2point[j][1];
      int point2 = seg2point[j][2];
      // Checking if the segment lies along one of the X cut planes or super-imposed exterior grid
      if (xcoord[point1-1]==X[i] && xcoord[point2-1]==X[i])
      {
        //checking that the segment in question is not equivalent to the entire grid segment. If it is equivalent, nothing is done, and the segment is preserved.
        if(((ycoord[point1-1]!=ymin)&&(ycoord[point2-1]!=ymax)) || ((ycoord[point2-1]!=ymin)&&(ycoord[point1-1]!=ymax)))
        {
          potentials.push_back(seg2point[j]);
        }				
      }
    }
    if (potentials.size() > 1)
    {
      vector <vector<double> > sorted = SortSegments(potentials,ycoord);
      int numsortedseg = sorted.size()/2;
      for (int j = 0; j < numsortedseg; ++j)
      {
        if (j == 0)
        {
          if (sorted[j][0] == ymin)
          {
            newseg.push_back(seg2point.size()+1);
            newseg.push_back(sorted[j+1][1]);
            newseg.push_back(sorted[j+2][2]);
            if (numbdmseg > 0)
            {
              newseg.push_back(-1);
            }
            seg2point.push_back(newseg);
            newseg.clear();
          }
          else
          {
            xcoord.push_back(X[i]);
            ycoord.push_back(ymin);
            if (numatt > 0)
            {
              for (int a = 0; a < numatt; ++a)
              {
                newattrow.push_back(0);
              }
              ptatt.push_back(newattrow);
              newattrow.clear();
            }
            if (numbdm > 0)
            {
              ptbdm.push_back(-1);
            }
            int MinPointNumber = xcoord.size();
            newseg.push_back(seg2point.size()+1);
            newseg.push_back(MinPointNumber);
            newseg.push_back(sorted[j][1]);
            if (numbdmseg > 0)
            {
              newseg.push_back(-1);
            }
            seg2point.push_back(newseg);
            newseg.clear();
          }
        }
        else if (j == numsortedseg-1)
        {
          if (sorted[sorted.size()-1][1] != ymax)
          {
            xcoord.push_back(X[i]);
            ycoord.push_back(ymax);
            if (numatt > 0)
            {
              for (int a = 0; a < numatt; ++a)
              {
                newattrow.push_back(0);
              }
              ptatt.push_back(newattrow);
              newattrow.clear();
            }
            if (numbdm > 0)
            {
              ptbdm.push_back(-1);
            }
            int MaxPointNumber = xcoord.size();
            newseg.push_back(seg2point.size()+1);
            newseg.push_back(sorted[sorted.size()-1][1]);
            newseg.push_back(MaxPointNumber);
            if (numbdmseg > 0)
            {
              newseg.push_back(-1);
            }
            seg2point.push_back(newseg);
            newseg.clear();
          }
        }
        else
        {
          int count = 2;
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(sorted[j+count][1]);
          newseg.push_back(sorted[j+count+1][1]);
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
          count = count + 1;
        }
      }
    }
    else if (potentials.size() == 1)
    {
      int pt1 = potentials[0][1];
      int pt2 = potentials[0][2];
      double y1 = ycoord[pt1-1];
      double y2 = ycoord[pt2-1];
      if (y1 < y2)
      {
        if (y1 != ymin)
        {
          xcoord.push_back(X[i]);
          ycoord.push_back(ymin);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(xcoord.size());
          newseg.push_back(pt1);
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
        if (y2 != ymax)
        {
          xcoord.push_back(X[i]);
          ycoord.push_back(ymax);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(pt2);
          newseg.push_back(xcoord.size());
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
      }
      else
      {
        if (y2 != ymin)
        {
          xcoord.push_back(X[i]);
          ycoord.push_back(ymin);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(xcoord.size());
          newseg.push_back(pt2);
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
        if (y1 != ymax)
        {
          xcoord.push_back(X[i]);
          ycoord.push_back(ymax);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(pt1);
          newseg.push_back(xcoord.size());
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
      }
    }
    else
    {
      xcoord.push_back(X[i]);
      ycoord.push_back(ymin);
      xcoord.push_back(X[i]);
      ycoord.push_back(ymax);
      newseg.push_back(seg2point.size()+1);
      newseg.push_back(xcoord.size()-1);
      newseg.push_back(xcoord.size());
      if (numbdmseg > 0)
      {
        newseg.push_back(-1);
      }
      seg2point.push_back(newseg);
      newseg.clear();
    }
    potentials.clear();
  }
  
  for (int i = 0; i < Y.size(); ++i)
  {
    for (int j = 0; j < seg2point.size(); ++j)
    {
      int point1 = seg2point[j][1];
      int point2 = seg2point[j][2];
      // Checking if the segment lies along one of the X cut planes or super-imposed exterior grid
      if (ycoord[point1-1]==Y[i] && ycoord[point2-1]==Y[i])
      {
        //checking that the segment in question is not equivalent to the entire grid segment. If it is equivalent, nothing is done, and the segment is preserved.
        if(((xcoord[point1-1]!=xmin)&&(xcoord[point2-1]!=xmax)) || ((xcoord[point2-1]!=xmin)&&(xcoord[point1-1]!=xmax)))
        {
          potentials.push_back(seg2point[j]);
        }				
      }
    }
    if (potentials.size() > 1)
    {
      vector <vector<double> > sorted = SortSegments(potentials,xcoord);
      int numsortedseg = sorted.size()/2;
      for (int j = 0; j < numsortedseg; ++j)
      {
        if (j == 0)
        {
          if (sorted[j][0] == xmin)
          {
            newseg.push_back(seg2point.size()+1);
            newseg.push_back(sorted[j+1][1]);
            newseg.push_back(sorted[j+2][2]);
            if (numbdmseg > 0)
            {
              newseg.push_back(-1);
            }
            seg2point.push_back(newseg);
            newseg.clear();
          }
          else
          {
            xcoord.push_back(xmin);
            ycoord.push_back(Y[i]);
            if (numatt > 0)
            {
              for (int a = 0; a < numatt; ++a)
              {
                newattrow.push_back(0);
              }
              ptatt.push_back(newattrow);
              newattrow.clear();
            }
            if (numbdm > 0)
            {
              ptbdm.push_back(-1);
            }
            int MinPointNumber = xcoord.size();
            newseg.push_back(seg2point.size()+1);
            newseg.push_back(MinPointNumber);
            newseg.push_back(sorted[j][1]);
            if (numbdmseg > 0)
            {
              newseg.push_back(-1);
            }
            seg2point.push_back(newseg);
            newseg.clear();
          }
        }
        else if (j == numsortedseg-1)
        {
          if (sorted[sorted.size()-1][1] != xmax)
          {
            xcoord.push_back(xmax);
            ycoord.push_back(Y[i]);
            if (numatt > 0)
            {
              for (int a = 0; a < numatt; ++a)
              {
                newattrow.push_back(0);
              }
              ptatt.push_back(newattrow);
              newattrow.clear();
            }
            if (numbdm > 0)
            {
              ptbdm.push_back(-1);
            }
            int MaxPointNumber = xcoord.size();
            newseg.push_back(seg2point.size()+1);
            newseg.push_back(sorted[sorted.size()-1][1]);
            newseg.push_back(MaxPointNumber);
            if (numbdmseg > 0)
            {
              newseg.push_back(-1);
            }
            seg2point.push_back(newseg);
            newseg.clear();
          }
        }
        else
        {
          int count = 2;
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(sorted[j+count][1]);
          newseg.push_back(sorted[j+count+1][1]);
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
          count = count + 1;
        }
      }
    }
    else if (potentials.size() == 1)
    {
      int pt1 = potentials[0][1];
      int pt2 = potentials[0][2];
      double x1 = xcoord[pt1-1];
      double x2 = xcoord[pt2-1];
      if (x1 < x2)
      {
        if (x1 != ymin)
        {
          xcoord.push_back(xmin);
          ycoord.push_back(Y[i]);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(xcoord.size());
          newseg.push_back(pt1);
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
        if (x2 != xmax)
        {
          xcoord.push_back(xmax);
          ycoord.push_back(Y[i]);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(pt2);
          newseg.push_back(xcoord.size());
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
      }
      else
      {
        if (x2 != xmin)
        {
          xcoord.push_back(xmin);
          ycoord.push_back(Y[i]);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(xcoord.size());
          newseg.push_back(pt2);
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
        if (x1 != xmax)
        {
          xcoord.push_back(xmax);
          ycoord.push_back(Y[i]);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              newattrow.push_back(0);
            }
            ptatt.push_back(newattrow);
            newattrow.clear();
          }
          if (numbdm > 0)
          {
            ptbdm.push_back(-1);
          }
          newseg.push_back(seg2point.size()+1);
          newseg.push_back(pt1);
          newseg.push_back(xcoord.size());
          if (numbdmseg > 0)
          {
            newseg.push_back(-1);
          }
          seg2point.push_back(newseg);
          newseg.clear();
        }
      }
    }
    else
    {
      xcoord.push_back(xmin);
      ycoord.push_back(Y[i]);
      xcoord.push_back(xmax);
      ycoord.push_back(Y[i]);
      newseg.push_back(seg2point.size()+1);
      newseg.push_back(xcoord.size()-1);
      newseg.push_back(xcoord.size());
      if (numbdmseg > 0)
      {
        newseg.push_back(-1);
      }
      seg2point.push_back(newseg);
      newseg.clear();
    }
    potentials.clear();
  }
}

int Tricall_Pre_Mesh_Gen::WhichSide(double Xmin, double Xmax, double Ymin, double Ymax, double x, double y)
{
/*Determines which side of a subdomain the two points of an edge will cut through. Returns a side number: 1 corresponds to bottom, 2 corresponds to right, 3 corresponds to top, and 4 corresponds to left. This function is only used in the case that one point is inside the cell and one point is outside it. Thus, it is only going to determine one side, because two cuts aren't possible. Xmin, Xmax, Ymin, and Ymax are the boundaries forming the subdomain, and x and y are the coordinates of the point.*/
  int side;
  if (x>=Xmin && x<=Xmax && y<=Ymin)
  {side = 1;}
  else if (y>=Ymin && y<=Ymax && x>=Xmax)
  {side = 2;}
  else if (x>=Xmin && x<=Xmax && y>=Ymax)
  {side = 3;}
  else if (y>=Ymin && y<=Ymax && x<=Xmin)
  {side = 4;}
  else
  {
    if (x>=Xmax && y<=Ymin)
    {
      if(abs(x-Xmax) <= abs(y-Ymin))
      {side = 1;}
      else
      {side = 2;}
    }
    else if (x<=Xmin && y<=Ymin)
    {
      if(abs(x-Xmin) <= abs(y-Ymin))
      {side = 1;}
      else
      {side = 4;}
    }
    else if (x>=Xmax && y>=Ymax)
    {
      if(abs(x-Xmax) <= abs(y-Ymax))
      {side = 3;}
      else
      {side = 2;}
    }
    else if (x<=Xmin && y>=Ymax)
    {
      if(abs(x-Xmin) <= abs(y-Ymax))
      {side = 3;}
      else
      {side = 4;}
    }
  }
  return side;
}

vector <int> Tricall_Pre_Mesh_Gen::WhichTwo(double Xmin,double Xmax,double Ymin,double Ymax,double x1,double x2,double y1,double y2)
{
//In the case where both points are outside, the edge could cut two sides or no sides. Input: the user-defined X and Y sections that serve to cut the input poly file; the two points that create the edge that is sectioned; i and j are the indices that define where in X and Y the function is looking at.

//This flag keeps track of the number of sides that intersect the edge outside of the cell. If the flag is equal to four, that means all intersection points are outside the cell, and therefore the edge does not cut any sides of the cell and is entirely outside.
int flag = 0;
//There are four possible "virtual" intersections with the cell. Because both points are outside, we need to figure out if they cut through the cell or not.

//This first case is done to determine and make sure that the edge does not intersect the grid in any way. Because there are four possible sides the grid could intersect, we must disprove all four sides.
  vector <int> sides;
  double x3; double y3;

  y3 = Yintersect(x1,y1,x2,y2,Xmin);
  if (y3>Ymax || y3<Ymin)
  {
    flag = flag + 1;
  }
  y3 = Yintersect(x1,y1,x2,y2,Xmax);
  if (y3>Ymax || y3<Ymin)
  {
    flag = flag + 1;
  }
  x3 = Xintersect(x1,y1,x2,y2,Ymin);
  if (x3>Xmax || x3<Xmin)
  {
    flag = flag + 1;
  }
  x3 = Xintersect(x1,y1,x2,y2,Ymax);
  if (x3>Xmax || x3<Xmin)
  {
    flag = flag + 1;
  }
      
  // if all sides have no true intersection point, then WhichTwo will return sides = [0 0];
  if (flag == 4)
  {
    sides.push_back(0);
    sides.push_back(0);
  }
  else
  {
    int side1 = WhichSide(Xmin,Xmax,Ymin,Ymax,x1,y1);
    sides.push_back(side1);
    int side2 = WhichSide(Xmin,Xmax,Ymin,Ymax,x2,y2);
    sides.push_back(side2);
  }
  return sides;
}

// This function finds out if a point is actually inside a cell or not. The point is entered in as a 1x2 row vector. If the point is inside, it sets a flag equal to 1. If not, flag = 0;
int Tricall_Pre_Mesh_Gen::IsInside(double x, double y, double Xmin, double Xmax, double Ymin, double Ymax)
{
  int flag;
  if (x>=Xmin && y>=Ymin && x<=Xmax && y<=Ymax)
  {
    flag = 1;
  }
  else
  {
    flag = 0;
  }
  return flag;
} 

void Tricall_Pre_Mesh_Gen::SegmentSplit(string & filename, double & Xmin, double & Xmax, double & Ymin, double & Ymax, vector <vector <double> > & seg2point, vector <double> & xcoord, vector <double> & ycoord, vector <vector <double> > & ptatt, vector <double> & ptbdm, int & numatt, int & numbdm, int & numbdmseg, vector <double> & newxcoord, vector <double> & newycoord, vector <vector <double> > & newptatt, vector <double> & newptbdm, vector <vector <double> > & SubDomSegments)
{
  int numpts = xcoord.size();
  int numseg = seg2point.size();
  vector <double> attrow;
  // a counter to keep track of the segments created
  int segcount = 1;
  for (int i = 0; i < numseg; ++i)
  {
    // a counter to keep track of the number of points created
    int countpts = 0;
    vector <double> CurrentSeg = seg2point[i];
    vector <double> newseg;
    int point1 = seg2point[i][1];
    int point2 = seg2point[i][2];
    // The x and y coordinates for the current segment
    double x1 = xcoord[point1-1];
    double x2 = xcoord[point2-1];
    double y1 = ycoord[point1-1];
    double y2 = ycoord[point2-1];
    // Checking whether the points are inside or outside of the current subdomain. If the flag is 0, the point is outside the subdomain, and if the flag is 1, the point is inside.
    int flag1 = IsInside(x1,y1,Xmin,Xmax,Ymin,Ymax);
    int flag2 = IsInside(x2,y2,Xmin,Xmax,Ymin,Ymax);
    // If the segment exists entirely within the subdomain, it is preserved
    if (flag1==1 && flag2==1)
    {
      newxcoord.push_back(x1); newycoord.push_back(y1);
      newxcoord.push_back(x2); newycoord.push_back(y2);
      if (numatt > 0)
      {
        for (int a = 0; a < numatt; ++a)
        {
          attrow.push_back(0);
        }
        newptatt.push_back(attrow);
        newptatt.push_back(attrow);
        attrow.clear();
      }
      if (numbdm > 0)
      {
        newptbdm.push_back(0);
        newptbdm.push_back(0);
      }
      newseg.push_back(segcount);
      newseg.push_back(xcoord.size());
      newseg.push_back(xcoord.size()-1);
      if (numbdmseg > 0)
      {
        newseg.push_back(CurrentSeg[3]);
      }
      SubDomSegments.push_back(newseg);
      newseg.clear();
      segcount = segcount + 1;
    }
    // otherwise, we must create points to shorten the segment so it exists within the subdomain
    else 
    {
      if (flag1==1 && flag2==0)
      {
        //Use WhichSide on the second point (since it is outside) to determine which side is being cut.
        int side = WhichSide(Xmin,Xmax,Ymin,Ymax,x2,y2);
        newxcoord.push_back(x1); newycoord.push_back(y1);
        if (numatt > 0)
        {
          for (int j = 0; j < numatt; ++j)
          {
            attrow.push_back(0);
          }
          newptatt.push_back(attrow);
        }
        if (numbdm > 0)
        {
          newptbdm.push_back(-1);
        }
        // the segment cuts the bottom subdomain boundary
        if (side == 1)
        {
          double x3min = Xintersect(x1,y1,x2,y2,Ymin);
          newxcoord.push_back(x3min);
          newycoord.push_back(Ymin);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        // the segment cuts the right subdomain boundary
        else if (side == 2)
        {
          double y3max = Yintersect(x1,y1,x2,y2,Xmax);
          newxcoord.push_back(Xmax);
          newycoord.push_back(y3max);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        // the segment cuts the top subdomain boundary
        else if (side == 3)
        {
          double x3max = Xintersect(x1,y1,x2,y2,Ymax);
          newxcoord.push_back(x3max);
          newycoord.push_back(Ymax);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        else
        {
          double y3min = Yintersect(x1,y1,x2,y2,Xmin);
          newxcoord.push_back(Xmin);
          newycoord.push_back(y3min);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        newseg.push_back(segcount);
        newseg.push_back(newxcoord.size());
        newseg.push_back(newxcoord.size()-1);
        if (numbdmseg > 0)
        {
          newseg.push_back(CurrentSeg[3]);
        }
        SubDomSegments.push_back(newseg);
        newseg.clear();
        segcount = segcount + 1;
      }
      else if (flag1==0 && flag2==1)
      {
        newxcoord.push_back(x2); newycoord.push_back(y2);
        if (numatt > 0)
        {
          for (int j = 0; j < numatt; ++j)
          {
            attrow.push_back(0);
          }
          newptatt.push_back(attrow);
        }
        if (numbdm > 0)
        {
          newptbdm.push_back(-1);
        }
        //Use WhichSide on the first point (since it is outside) to determine which side is being cut.
        int side = WhichSide(Xmin,Xmax,Ymin,Ymax,x1,y1);
        // the segment cuts the bottom subdomain boundary
        if (side == 1)
        {
          double x3min = Xintersect(x1,y1,x2,y2,Ymin);
          newxcoord.push_back(x3min);
          newycoord.push_back(Ymin);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        // the segment cuts the right subdomain boundary
        else if (side == 2)
        {
          double y3max = Yintersect(x1,y1,x2,y2,Xmax);
          newxcoord.push_back(Xmax);
          newycoord.push_back(y3max);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        // the segment cuts the top subdomain boundary
        else if (side == 3)
        {
          double x3max = Xintersect(x1,y1,x2,y2,Ymax);
          newxcoord.push_back(x3max);
          newycoord.push_back(Ymax);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        else
        {
          double y3min = Yintersect(x1,y1,x2,y2,Xmin);
          newxcoord.push_back(Xmin);
          newycoord.push_back(y3min);
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
          }
        }
        newseg.push_back(segcount);
        newseg.push_back(newxcoord.size());
        newseg.push_back(newxcoord.size()-1);
        if (numbdmseg > 0)
        {
          newseg.push_back(CurrentSeg[3]);
        }
        SubDomSegments.push_back(newseg);
        newseg.clear();
        segcount = segcount + 1;
      }
      // Both points are outside
      else
      {
        //Use WhichTwo to check which two boundaries the segment cuts through
        vector <int> sides;
        sides = WhichTwo(Xmin,Xmax,Ymin,Ymax,x1,x2,y1,y2);
        int side1 = sides[0]; int side2 = sides[1];
        if (side1!=0 && side2!=0)
        {
          if ((side1==1 && side2==2) || (side1==2 && side2==1))
          {
            double y3 = Ymin;
            double x3 = Xintersect(x1,y1,x2,y2,y3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
            
            x3 = Xmax;
            y3 = Yintersect(x1,y1,x2,y2,x3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
          }
          else if ((side1==1 && side2==3) || (side1==3 && side2==1))
          {
            double y3 = Ymin;
            double x3 = Xintersect(x1,y1,x2,y2,y3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
            
            y3 = Ymax;
            x3 = Xintersect(x1,y1,x2,y2,y3);
            newxcoord.push_back(x3); newycoord.push_back(y3);	
          }
          else if ((side1==1 && side2==4) || (side1==4 && side2==1))
          {
            double y3 = Ymin;
            double x3 = Xintersect(x1,y1,x2,y2,y3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
            
            x3 = Xmin;
            y3 = Yintersect(x1,y1,x2,y2,x3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
          }
          else if ((side1==2 && side2==3) || (side1==3 && side2==2))
          {
            double y3 = Ymax;
            double x3 = Xintersect(x1,y1,x2,y2,y3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
            
            x3 = Xmax;
            y3 = Yintersect(x1,y1,x2,y2,x3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
          }
          else if ((side1==2 && side2==4) || (side1==4 && side2==2))
          {
            double x3 = Xmin;
            double y3 = Yintersect(x1,y1,x2,y2,x3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
            
            x3 = Xmax;
            y3 = Yintersect(x1,y1,x2,y2,x3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
          }
          else if ((side1==3 && side2==4) || (side1==4 && side2==3))
          {
            double x3 = Xmin;
            double y3 = Yintersect(x1,y1,x2,y2,x3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
            
            y3 = Ymax;
            x3 = Xintersect(x1,y1,x2,y2,y3);
            newxcoord.push_back(x3); newycoord.push_back(y3);
          }
          // adding in attributes and boundary markers for the two new points
          if (numatt > 0)
          {
            for (int a = 0; a < numatt; ++a)
            {
              attrow.push_back(0);
            }
            newptatt.push_back(attrow);
            newptatt.push_back(attrow);
            attrow.clear();
          }
          if (numbdm > 0)
          {
            newptbdm.push_back(0);
            newptbdm.push_back(0);
          }
          newseg.push_back(segcount);
          newseg.push_back(newxcoord.size());
          newseg.push_back(newxcoord.size()-1);
          if (numbdmseg > 0)
          {
            newseg.push_back(CurrentSeg[3]);
          }
          SubDomSegments.push_back(newseg);
          newseg.clear();
          segcount = segcount + 1;
        }
      }
    }
  }
  // Making sure no points are repeated
  // Keeps track of which points were erased
  map <int,int> erased;
  // Keeping track of how many points were erased
  int erasecount = 0;
  // Two vectors that preserve the original x and y coordinates before making them unique
  vector <double> holdx = newxcoord; vector <double> holdy = newycoord;
  for (int i = 0; i < newxcoord.size(); ++ i)
  {
    double x = newxcoord[i]; double y = newycoord[i];
    for (int j = i+1; j < newxcoord.size(); ++j)
    {
      if (x==newxcoord[j] && y == newycoord[j])
      {
        newxcoord.erase(newxcoord.begin()+j);
        newycoord.erase(newycoord.begin()+j);
      }
    }
  }

  // Renumbering the points in SubDomSegments so that they correspond to SubDomSegments
  for (int i = 0; i < SubDomSegments.size(); ++i)
  {
    int pt1 = SubDomSegments[i][1];
    int pt2 = SubDomSegments[i][2];
    double x = holdx[pt1-1]; double y = holdy[pt1-1];
    for (int j = 0; j < newxcoord.size(); ++j)
    {
      if (x == newxcoord[j] && y == newycoord[j])
      {
        SubDomSegments[i][1] = j+1;
      }
    }

    x = holdx[pt2-1]; y = holdy[pt2-1];
    for (int j = 0; j < newxcoord.size(); ++j)
    {
      if (x == newxcoord[j] && y == newycoord[j])
      {
        SubDomSegments[i][2] = j+1;
      }
    }
  }

}

int Tricall_Pre_Mesh_Gen::Neigh(int NeighX, int NeighY, map <string, int> SubCoordToInt)
{
  stringstream coordconv;
  coordconv << NeighX << ',' << NeighY;
  string NeighCoord = coordconv.str();
  int NeighID = SubCoordToInt[NeighCoord];
  
  return NeighID;
}


void Tricall_Pre_Mesh_Gen::Centroid(vector <vector<double> > & polygonpts, double & Cx, double & Cy)
{
  // We need to make each point unique 
  for (int i = 0; i < polygonpts.size(); ++i)
  {
    double x = polygonpts[i][0]; double y = polygonpts[i][1];
    for (int j = i+1; j < polygonpts.size()-1; ++j)
    {
      if (x==polygonpts[j][0] && y==polygonpts[j][1])
      {
        polygonpts.erase(polygonpts.begin()+j);
      }
    }
  }

  // We need to sort the vertices in a counter-clockwise order
  vector <vector <double> > sortedpts;
  vector <double> row;
  vector < vector <double> > polarAngle;
  // Beginning polar coord sort experiment
  // First, we need to find the lowest y-coordinate to use as our reference point
  sort(polygonpts.begin(),polygonpts.end(),sortFunc);
  // The reference coordinates
  double baseX = polygonpts[0][0]; double baseY = polygonpts[0][1];
  // Finding the polar angle for each coordinate using the above coordinates as the reference
  for (int i = 1; i < polygonpts.size(); ++i)
  {
    double x = polygonpts[i][0]-baseX; double y = polygonpts[i][1]-baseY;
    double phi = atan2(y,x);
    // Storing the original index and the polar angle
    row.push_back(i); row.push_back(phi);
    polarAngle.push_back(row); row.clear();
  }
  sort(polarAngle.begin(),polarAngle.end(),sortFunc);
  
  // Putting the reference point into sortedpts
  row.push_back(baseX); row.push_back(baseY);
  sortedpts.push_back(row); row.clear();
  
  for (int i = 0; i < polarAngle.size(); ++i)
  {
    int index = polarAngle[i][0];
    row.push_back(polygonpts[index][0]); row.push_back(polygonpts[index][1]);
    sortedpts.push_back(row); row.clear();
  }
  polygonpts = sortedpts;
  polygonpts.push_back(polygonpts[0]);

  // The area of the polygon
  double Area;
  // A summation used to compute the area
  double AreaSum = 0;
  // The number of vertices in the polygon
  int n = polygonpts.size();
  
  vector <double> x;
  vector <double> y; 
  for (int i = 0; i < polygonpts.size(); ++i)
  {
    x.push_back(polygonpts[i][0]);
    y.push_back(polygonpts[i][1]);
  }
  // Calculating the area
  for (int i = 0; i <= n-1; ++i)
  {
    AreaSum = AreaSum + (x[i]*y[i+1] - x[i+1]*y[i]);
  }
  Area = AreaSum/2;
  
  double Xsum = 0;
  double Ysum = 0;
  for (int i = 0; i <= n-1; ++i)
  {
    Xsum = Xsum + (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i]);
    Ysum = Ysum + (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i]);
  }
  Cx = Xsum/(6*Area);
  Cy = Ysum/(6*Area);
  
}

vector <vector <double> > Tricall_Pre_Mesh_Gen::RegionalAtt(string fileprefix,double Xmin,double Xmax,double Ymin,double Ymax)
{
  // the number of triangles
  int ntri;
  // the regional attributes propagated in each triangle
  vector <double> AttID;
  // the edge numbers that create each triangle
  vector < vector <int> > tri2edge;
  // the two triangle that create each edge
  vector < vector <int> > edge2tri;
  // the points creating each edge
  vector < vector <int> > edge2pt;
  // the points corresponding to each triangle
  vector < vector <int> > tri2pt;
  // the coordinates corresponding to each point
  vector < vector <double> > pt2coord;
  // which output we are looking at
  int index = 1;
  // compiling information from the constrained Delaunay Triangulation
  load_geom(ntri,AttID,tri2edge,edge2tri,edge2pt,tri2pt,pt2coord,fileprefix,index);

  vector <double> xcoord;
  for (int i = 0; i < pt2coord.size(); ++i)
  {
    xcoord.push_back(pt2coord[i][0]);
  }
  vector <double> ycoord;
  for (int i = 0; i < pt2coord.size(); ++i)
  {
    ycoord.push_back(pt2coord[i][1]);
  }
  // A vector storing the points of the polygon formed by the triangle being portioned
  vector < vector <double> > polygonpts;
  // A vector that stores attribute information
  vector <vector <double> > RegionalAtt;
  vector <double> regionrow;
  int side;
  for (int a = 0; a < ntri; ++a)
  {
    // the regional attribute for the current triangle
    double triID = AttID[a];
    // this triangles edges
    int edge1 = tri2edge[a][0];
    int edge2 = tri2edge[a][1];
    int edge3 = tri2edge[a][2];
    vector <int> edges;
    edges.push_back(edge1); edges.push_back(edge2); edges.push_back(edge3);
    // looping through the edges of the current triangle 
    for (int k = 0; k < edges.size(); ++k)
    {
      int CurrentEdge = edges[k];
      int endpt1 = edge2pt[CurrentEdge-1][0];
      int endpt2 = edge2pt[CurrentEdge-1][1];
      double x1 = xcoord[endpt1-1]; double y1 = ycoord[endpt1-1];
      double x2 = xcoord[endpt2-1]; double y2 = ycoord[endpt2-1];
      // Checking to see if the two points in the edge are in the same cell or not
      int flag1 = IsInside(x1,y1,Xmin,Xmax,Ymin,Ymax);
      int flag2 = IsInside(x2,y2,Xmin,Xmax,Ymin,Ymax);
      // Temp vectors to store info that will later go into polygonpts
      vector <double> row1; vector <double> row2;
      // If both points of the edge are in the same cell
      if (flag1==1 && flag2==1)
      {
        row1.push_back(x1); row1.push_back(y1);
        row2.push_back(x2); row2.push_back(y2);
        polygonpts.push_back(row1); row1.clear();
        polygonpts.push_back(row2); row2.clear();
      }
      // If not, we must create a point to shorten the edge so it exists within the subdomain
      else
      {
        // Endpoint1 is in the subdomain, and Endpoint2 is outside
        if ((flag1==1 && flag2==0) || (flag1==0 && flag2==1))
        {
          if (flag1 == 1)
          {
            // WhichSide will tell us where the second point is, so that we know which side of the subdomain is being cut by the edge
            side = WhichSide(Xmin,Xmax,Ymin,Ymax,x2,y2);
            // Adding point1 into the polygon
            row1.push_back(x1); row1.push_back(y1);
            polygonpts.push_back(row1); row1.clear();
          }
          if (flag2 == 1)
          {
            // WhichSide will tell us where the first point is, so that we know which side of the subdomain is being cut by the edge
            side = WhichSide(Xmin,Xmax,Ymin,Ymax,x1,y1);
            // Adding point1 into the polygon
            row1.push_back(x2); row1.push_back(y2);
            polygonpts.push_back(row1); row1.clear();
          }
          
          // The edge cuts the bottom subdomain boundary
          if (side == 1)
          {
            double x3min = Xintersect(x1,y1,x2,y2,Ymin);
            row2.push_back(x3min); row2.push_back(Ymin);
            polygonpts.push_back(row2); row2.clear();
          }
          // The edge cuts the right subdomain boundary
          else if (side == 2)
          {
            double y3max = Yintersect(x1,y1,x2,y2,Xmax);
            row2.push_back(Xmax); row2.push_back(y3max);
            polygonpts.push_back(row2); row2.clear();
          }
          // The edge cuts the top subdomain boundary
          else if (side == 3)
          {
            double x3max = Xintersect(x1,y1,x2,y2,Ymax);
            row2.push_back(x3max); row2.push_back(Ymax);
            polygonpts.push_back(row2); row2.clear();
          }
          // The edge cuts the left subdomain boundary
          else
          {
            double y3min = Yintersect(x1,y1,x2,y2,Xmin);
            row2.push_back(Xmin); row2.push_back(y3min);
            polygonpts.push_back(row2); row2.clear();
          }
        }
        // If both points of the edge are outside, we need to check if they cut either two or none of the subdomain boundaries
        else
        {
          // Use WhichTwo to check which boundaries (if any) are cut
          vector <int> sides;
          sides = WhichTwo(Xmin,Xmax,Ymin,Ymax,x1,x2,y1,y2);
          int side1 = sides[0]; int side2 = sides[1];			
          if (side1!=0 && side2!=0)
          {
            if ((side1==1 && side2==2) || (side1==2 && side2==1))
            {
              double y3 = Ymin;
              double x3 = Xintersect(x1,y1,x2,y2,y3);
              row1.push_back(x3); row1.push_back(y3);
              polygonpts.push_back(row1); row1.clear();
              
              x3 = Xmax;
              y3 = Yintersect(x1,y1,x2,y2,x3);
              row2.push_back(x3); row2.push_back(y3);
              polygonpts.push_back(row2); row2.clear();				
            }
            else if ((side1==1 && side2==3) || (side1==3 && side2==1))
            {
              double y3 = Ymin;
              double x3 = Xintersect(x1,y1,x2,y2,y3);
              row1.push_back(x3); row1.push_back(y3);
              polygonpts.push_back(row1); row1.clear();
              
              y3 = Ymax;
              x3 = Xintersect(x1,y1,x2,y2,y3);
              row2.push_back(x3); row2.push_back(y3);
              polygonpts.push_back(row2); row2.clear();
            }
            else if ((side1==1 && side2==4) || (side1==4 && side2==1))
            {
              double y3 = Ymin;
              double x3 = Xintersect(x1,y1,x2,y2,y3);
              row1.push_back(x3); row1.push_back(y3);
              polygonpts.push_back(row1); row1.clear();
              
              x3 = Xmin;
              y3 = Yintersect(x1,y1,x2,y2,x3);
              row2.push_back(x3); row2.push_back(y3);
              polygonpts.push_back(row2); row2.clear();
            }
            else if ((side1==2 && side2==3) || (side1==3 && side2==2))
            {
              double y3 = Ymax;
              double x3 = Xintersect(x1,y1,x2,y2,y3);
              row1.push_back(x3); row1.push_back(y3);
              polygonpts.push_back(row1); row1.clear();
              
              x3 = Xmax;
              y3 = Yintersect(x1,y1,x2,y2,x3);
              row2.push_back(x3); row2.push_back(y3);
              polygonpts.push_back(row2); row2.clear();
            }
            else if ((side1==2 && side2==4) || (side1==4 && side2==2))
            {
              double x3 = Xmin;
              double y3 = Yintersect(x1,y1,x2,y2,x3);
              row1.push_back(x3); row1.push_back(y3);
              polygonpts.push_back(row1); row1.clear();
              
              x3 = Xmax;
              y3 = Yintersect(x1,y1,x2,y2,x3);
              row2.push_back(x3); row2.push_back(y3);
              polygonpts.push_back(row2); row2.clear();
            }
            else if ((side1==3 && side2==4) || (side1==4 && side2==3))
            {
              double x3 = Xmin;
              double y3 = Yintersect(x1,y1,x2,y2,x3);
              row1.push_back(x3); row1.push_back(y3);
              polygonpts.push_back(row1); row1.clear();
              
              y3 = Ymax;
              x3 = Xintersect(x1,y1,x2,y2,y3);
              row2.push_back(x3); row2.push_back(y3);
              polygonpts.push_back(row2); row2.clear();
            }
          }
        }
      }	
    }
    // Finding the Centroid of the newly created polygon.
    if (polygonpts.size() > 0 )
    {
      double Cx; double Cy;
      Centroid(polygonpts,Cx,Cy);
      regionrow.push_back(Cx); 
      regionrow.push_back(Cy); 
      regionrow.push_back(triID);
      RegionalAtt.push_back(regionrow);
      polygonpts.clear();
      regionrow.clear();
    }
  }
  return RegionalAtt;
}
 

struct main_wf
{
  int    m_argc;
  char** m_argv; // Modification of argv is not permitted
  MPI_Comm m_comm;

  main_wf(int argc, char** argv, MPI_Comm comm)
    : m_argc(argc), m_argv(argv), m_comm(comm)
  {}

  void operator()(void)
  {
  
  Tricall_Utility tu;
  //tu.trianglePreload();
  Tricall_Pre_Mesh_Gen tricall;
  tricall.generateMesh();
  
   
  }
};

int main(int argc, char** argv)
{
  int mpi_init;
  MPI_Init(&argc, &argv);
  MPI_Initialized(&mpi_init);
  if(!mpi_init)
  {
    std::cerr << "MPI initialization failed." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  stapl::option opt(argc, argv);
  opt = opt & stapl::option("MPI_Comm", MPI_Comm(MPI_COMM_WORLD));
  stapl::initialize(opt);
  stapl::execute(::main_wf(argc, argv, MPI_COMM_WORLD),
                 opt.get<unsigned int>("STAPL_MAIN_LEVELS", 1));
  stapl::finalize();

  MPI_Finalize();

  return EXIT_SUCCESS;
}




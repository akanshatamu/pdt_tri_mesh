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
#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

////////////////////////////////////////////////////////////////////////////////
///
/// @brief connect to an external object file and declare the methods to be 
///  called in this file.
///
////////////////////////////////////////////////////////////////////////////////
extern "C" {
    #include "triangle.h"
	void triangulate(char *, struct triangulateio *, struct triangulateio *,
                 struct triangulateio *,char *, int);
	//void writepoly(struct mesh *m, struct behavior *b, char *polyfilename,
     //          REAL *holelist, int holes, REAL *regionlist, int regions,
      //         int argc, char **argv);
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief .
///
/// This class is used to 1> Call triangle mesh generator library 2> export 
/// mesh detail to an external file.
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "Data_Transfer_Object.h"

using std::cout;
using std::ofstream;
using std::ios;
using std::vector;
using std::string;
using std::ostringstream;

class Tricall_Utility
{
private:
	//

public:
	Tricall_Utility(void);
	~Tricall_Utility(void);

	void report(struct triangulateio *io, int markers, int reporttriangles, int reportneighbors, int reportsegments,
            int reportedges, int reportnorms);
	void exportnodes(struct triangulateio *io,int markers,char *filename);
	void exportedges(struct triangulateio *io,int markers,char *filename);
	void exportelements(struct triangulateio *io,char *filename);
	void printHowdy();
	void trianglePreload();
	void generateTriangleMesh(Data_Transfer_Object&);
	void createPolyFiles(vector < vector < vector <double> > >,int,int,int,vector< string > );
};


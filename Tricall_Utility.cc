#include "Tricall_Utility.h"


/**
'z' switch so that points are numbered from 0
'N' don't need the boundary markers

 **/


Tricall_Utility::Tricall_Utility(void)
{
}


Tricall_Utility::~Tricall_Utility(void)
{
}


void Tricall_Utility::printHowdy()
{
	std::cout << "Howdy From Tricall_Utility" ;
}

//This is currently used
void Tricall_Utility::generateTriangleMesh(Data_Transfer_Object &dto)
{
	std::cout << "Howdy From generateTriangleMesh\n" ;
	
	vector<double> pointlist=dto.pointlist;
	vector<double> pointattributelist=dto.pointattributelist;
	vector<int> pointmarkerlist=dto.pointmarkerlist;
	vector<double> regionlist=dto.regionlist;
	vector<int> segmentlist=dto.segmentlist;
	vector<int> segmentmarkerlist=dto.segmentmarkerlist;
	double maximumTriangleArea=dto.maximumTriangleArea;
	
	Tricall_Utility tri;
	
	struct triangulateio in, mid, out, vorout;

  /* Define input points. */

  in.numberofpoints = dto.numberofpoints;
  in.numberofpointattributes = dto.numberofpointattributes;
  
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));  
  for(int iTer=0;iTer < (in.numberofpoints * 2);iTer++){
	in.pointlist[iTer]=pointlist[iTer];
  }
  
  in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                          in.numberofpointattributes *
                                          sizeof(REAL));
  for(int iTer=0;iTer < (in.numberofpoints*in.numberofpointattributes);iTer++){
	in.pointattributelist[iTer]=pointattributelist[iTer];
  }									  
										  
  in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
  for(int iTer=0;iTer < (in.numberofpoints);iTer++){
	in.pointmarkerlist[iTer]=pointmarkerlist[iTer];
  }
  
  in.numberofsegments = dto.numberofsegments;
  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  for(int iTer=0;iTer < (in.numberofsegments * 2 );iTer++){
	in.segmentlist[iTer]=segmentlist[iTer];
  }
  in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
  for(int iTer=0;iTer < (in.numberofsegments);iTer++){
	in.segmentmarkerlist[iTer]=segmentmarkerlist[iTer];
  }
  
  in.numberofholes = dto.numberofholes;
  //in.holelist = (REAL *) malloc(in.numberofholes * 4 * sizeof(REAL));
  
  in.numberofregions = dto.numberofregions;  
  in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
  for(int iTer=0;iTer < (in.numberofregions * 4 );iTer++){
	in.regionlist[iTer]=regionlist[iTer];
  }

  printf("Input point set:\n\n");
  tri.report(&in, 1, 0, 0, 0, 0, 0);

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

  mid.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  mid.pointattributelist = (REAL *) NULL;
  mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  mid.triangleattributelist = (REAL *) NULL;
  mid.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  mid.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  mid.segmentmarkerlist = (int *) NULL;
  mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

  vorout.pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
  /* Needed only if -v switch used and number of attributes is not zero: */
  vorout.pointattributelist = (REAL *) NULL;
  vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
  vorout.normlist = (REAL *) NULL;         /* Needed only if -v switch used. */

  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), number everything from  */
  /*   zero (z), assign a regional attribute to each element (A), and  */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  char *testval="begin";
  triangulate("pcAevn", &in, &mid, &vorout, testval, 0);
  //triangulate("pcecAevn", &in, &mid, &vorout);

  //printf("Initial triangulation:\n\n");
  //tri.report(&mid, 1, 1, 1, 1, 1, 0);
  //printf("Initial Voronoi diagram:\n\n");
  //tri.report(&vorout, 0, 0, 0, 0, 1, 1);

  /* Attach area constraints to the triangles in preparation for */
  /*   refining the triangulation.                               */

  /* Needed only if -r and -a switches used: */
  mid.trianglearealist = (REAL *) malloc(mid.numberoftriangles * sizeof(REAL));
  mid.trianglearealist[0] = maximumTriangleArea;
  mid.trianglearealist[1] = maximumTriangleArea;

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out'.                                    */

  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of attributes is zero: */
  out.pointattributelist = (REAL *) NULL;
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  out.triangleattributelist = (REAL *) NULL;

  /* Refine the triangulation according to the attached */
  /*   triangle area constraints.                       */

  char *filename=dto.filenamestr;
  //triangulate("praBP", &mid, &out, (struct triangulateio *) NULL);
  triangulate("enpraDP", &mid, &out, (struct triangulateio *) NULL, filename, 1);

  //printf("Refined triangulation:\n\n");
  //tri.report(&out, 1, 1, 1, 0, 1, 0);
  
  tri.exportnodes(&out,1,filename);
  tri.exportedges(&out,1,filename);
  tri.exportelements(&out,filename);
  
  /* Free all allocated arrays, including those allocated by Triangle. */

  free(in.pointlist);
  free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.regionlist);
  free(mid.pointlist);
  free(mid.pointattributelist);
  free(mid.pointmarkerlist);
  free(mid.trianglelist);
  free(mid.triangleattributelist);
  free(mid.trianglearealist);
  free(mid.neighborlist);
  free(mid.segmentlist);
  free(mid.segmentmarkerlist);
  free(mid.edgelist);
  free(mid.edgemarkerlist);
  free(vorout.pointlist);
  free(vorout.pointattributelist);
  free(vorout.edgelist);
  free(vorout.normlist);
  free(out.pointlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);	
	
	
}

//Not currently using 
//This functionaly is used in tricall_pre_mesh_gen
void Tricall_Utility::createPolyFiles(vector < vector < vector <double> > > polycell, 
	int numcells,int numatt,int numbdm, vector <string> filenames )
{
	std::cout << "Howdy From createPolyFiles" ;
	
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
			fprintf(outfile, "%s\n", row.c_str());
			//newfile << row;
			//newfile<< endl;
		}
		//newfile.close();
		fclose(outfile);
	}
}

void Tricall_Utility::exportelements(struct triangulateio *io,char *filename)
{
	printf("Entered export elements\n");
	FILE *outfile;

	//Start concatenation
	//Concatenation two char* strings
	//To dynamically create the filename with ".node" appended 
	//to the end of the filename
	char *extension=".ele";
	int len1 = strlen(filename);
	int len2 = strlen(extension);
	char nodefilename[len1 + len2 + 1];
	strcpy(nodefilename,filename);
	strcat(nodefilename,extension);
	//End concatenation
	printf("-----------Writing %s.\n", nodefilename);

	//Open the file
	outfile = fopen(nodefilename, "w+");
	if (outfile == (FILE *) NULL) {
	printf("  Error:  Cannot create file %s.\n", nodefilename);
	}
	
	/* Number of triangles, vertices per triangle, attributes per triangle. */
	fprintf(outfile, "%d  %d  %d\n", io->numberoftriangles,
		  io->numberofcorners, io->numberofpointattributes);
	
	for (int i = 0; i < io->numberoftriangles; i++) {
		
		printf("Triangle %4d points:", i);
		fprintf(outfile, "%4ld    ", i);
		for (int j = 0; j < io->numberofcorners; j++) {
		  printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
		  fprintf(outfile, "%4d  ", io->trianglelist[i * io->numberofcorners + j]);
		}
		if (io->numberoftriangleattributes > 0) {
		  printf("   attributes");
		}
		for (int j = 0; j < io->numberoftriangleattributes; j++) {
		  printf("  %.6g", io->triangleattributelist[i *
										 io->numberoftriangleattributes + j]);
		  fprintf(outfile, "%.17g  ", io->triangleattributelist[i *
										 io->numberoftriangleattributes + j]);								 
		}
		printf("\n");
		fprintf(outfile, "\n");		
	}
	printf("\n");
	
}

void Tricall_Utility::exportedges(struct triangulateio *io, int markers, char *filename)
{
	printf("Entered export edges\n");
	FILE *outfile;

	//Start concatenation
	//Concatenation two char* strings
	//To dynamically create the filename with ".node" appended 
	//to the end of the filename
	char *extension=".edge";
	int len1 = strlen(filename);
	int len2 = strlen(extension);
	char nodefilename[len1 + len2 + 1];
	strcpy(nodefilename,filename);
	strcat(nodefilename,extension);
	//End concatenation
	printf("-----------Writing %s.\n", nodefilename);

	//Open the file
	outfile = fopen(nodefilename, "w+");
	if (outfile == (FILE *) NULL) {
	printf("  Error:  Cannot create file %s.\n", nodefilename);
	}
	
	/* Number of edges, number of boundary markers (zero or one). */
	fprintf(outfile, "%d  %d\n", io->numberofedges, io->bdr_markers);
	printf("Mid writing edge files\n");
	
	for (int i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
	  fprintf(outfile, "%4ld   ", i);
      for (int j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
		fprintf(outfile, "%d  ", io->edgelist[i * 2 + j]);
      }
      
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
		fprintf(outfile, "%d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
		fprintf(outfile, "\n");
      }
    }
    printf("\n");
	
}

void Tricall_Utility::exportnodes(struct triangulateio *io, int markers, char *filename)
{
  int i, j;

  printf("Entered export nodes\n");
  FILE *outfile;
  
  //Start concatenation
  //Concatenation two char* strings
  //To dynamically create the filename with ".node" appended 
  //to the end of the filename
  char *extension=".node";
  int len1 = strlen(filename);
  int len2 = strlen(extension);
  char nodefilename[len1 + len2 + 1];
  strcpy(nodefilename,filename);
  strcat(nodefilename,extension);
  //End concatenation
  printf("-----------Writing %s.\n", nodefilename);
  
  //Open the file
  outfile = fopen(nodefilename, "w+");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", nodefilename);
  }
  /* Number of vertices, number of dimensions, number of vertex attributes, */
  /*   and number of boundary markers (zero or one).                        */
  fprintf(outfile, "%d  %d  %d  %d\n", io->numberofpoints, io->mesh_dim,
          io->numberofpointattributes, io->bdr_markers);
  
  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
	
	/* Vertex number, x and y coordinates. */
    fprintf(outfile, "%4d    %.17g  %.17g", i, io->pointlist[i * 2 + 0],
              io->pointlist[i * 2 + 1]);
	
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
	  /* Write an attribute. */
      fprintf(outfile, "  %.17g", io->pointattributelist[i * io->numberofpointattributes + j]);		 
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
	  /* Write the boundary marker. */
      fprintf(outfile, "    %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
	  fprintf(outfile, "\n");
    }
  }
  printf("\n");
  fclose(outfile);
  printf("Exiting exportnodes\n");
}

void Tricall_Utility::report(struct triangulateio *io, int markers, int reporttriangles, int reportneighbors, int reportsegments,
            int reportedges, int reportnorms)
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

    
  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}

/**int main(int argc, char** argv)
{
   std::cout << "Howdy From main,tricall_utility" ;
   tricall_utility tri;
   tri.trianglePreload();
   
}**/

//With default values
//USed to test triangle code
//should be used as main with a different makefile
void Tricall_Utility::trianglePreload()
{

	std::cout << "trianglePreload" ;
	Tricall_Utility tri;
	
	struct triangulateio in, mid, out, vorout;

  /* Define input points. */

  in.numberofpoints = 6;
  in.numberofpointattributes = 0;
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
  in.pointlist[0] = 0.2;
  in.pointlist[1] = 0.0;
  in.pointlist[2] = 0.5;
  in.pointlist[3] = 0.0;
  in.pointlist[4] = 0.2;
  in.pointlist[5] = 0.0;
  in.pointlist[6] = 0.2;
  in.pointlist[7] = 0.5;
  in.pointlist[8] = 0.5;
  in.pointlist[9] = 0.0;
  in.pointlist[10] = 0.5;
  in.pointlist[11] = 0.5;
  
  in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                          in.numberofpointattributes *
                                          sizeof(REAL));
  in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
  in.pointmarkerlist[0] = 0;
  in.pointmarkerlist[1] = 1;
  in.pointmarkerlist[2] = 0;
  in.pointmarkerlist[3] = 2;
  in.pointmarkerlist[4] = 0;
  in.pointmarkerlist[5] = 0;

  
  in.numberofsegments = 0;
  in.numberofholes = 0;
  in.numberofregions = 0;
  in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));

  printf("Input point set:\n\n");
  tri.report(&in, 1, 0, 0, 0, 0, 0);

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

  mid.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  mid.pointattributelist = (REAL *) NULL;
  mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  mid.triangleattributelist = (REAL *) NULL;
  mid.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  mid.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  mid.segmentmarkerlist = (int *) NULL;
  mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

  vorout.pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
  /* Needed only if -v switch used and number of attributes is not zero: */
  vorout.pointattributelist = (REAL *) NULL;
  vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
  vorout.normlist = (REAL *) NULL;         /* Needed only if -v switch used. */

  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), number everything from  */
  /*   zero (z), assign a regional attribute to each element (A), and  */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  char *testval="begin";
  triangulate("pcAevn", &in, &mid, &vorout, testval, 1);
  //triangulate("pcecAevn", &in, &mid, &vorout);

  printf("Initial triangulation:\n\n");
  tri.report(&mid, 1, 1, 1, 1, 1, 0);
  printf("Initial Voronoi diagram:\n\n");
  tri.report(&vorout, 0, 0, 0, 0, 1, 1);

  /* Attach area constraints to the triangles in preparation for */
  /*   refining the triangulation.                               */

  /* Needed only if -r and -a switches used: */
  mid.trianglearealist = (REAL *) malloc(mid.numberoftriangles * sizeof(REAL));
  mid.trianglearealist[0] = 0.03;
  mid.trianglearealist[1] = 0.01;

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out'.                                    */

  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of attributes is zero: */
  out.pointattributelist = (REAL *) NULL;
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  out.triangleattributelist = (REAL *) NULL;

  /* Refine the triangulation according to the attached */
  /*   triangle area constraints.                       */

  char *testval1="refine";
  //triangulate("praBP", &mid, &out, (struct triangulateio *) NULL);
  triangulate("enpraDP", &mid, &out, (struct triangulateio *) NULL, testval1, 1);

  printf("Refined triangulation:\n\n");
  tri.report(&out, 1, 1, 1, 0, 1, 0);
  
  char *filenamestr="test2";
  //tri.exportnodes(&out,1,filenamestr);
  
  
  
  /* Free all allocated arrays, including those allocated by Triangle. */

  free(in.pointlist);
  free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.regionlist);
  free(mid.pointlist);
  free(mid.pointattributelist);
  free(mid.pointmarkerlist);
  free(mid.trianglelist);
  free(mid.triangleattributelist);
  free(mid.trianglearealist);
  free(mid.neighborlist);
  free(mid.segmentlist);
  free(mid.segmentmarkerlist);
  free(mid.edgelist);
  free(mid.edgemarkerlist);
  free(vorout.pointlist);
  free(vorout.pointattributelist);
  free(vorout.edgelist);
  free(vorout.normlist);
  free(out.pointlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);	
	
}







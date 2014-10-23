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


//This is currently used
void Tricall_Utility::generateMeshWithFixedBoundaries(Data_Transfer_Object &dto)
{
  std::cout << "Howdy From generateMeshWithFixedBoundaries for " << dto.filenamestr << "\n" ;
  vector <double> edgelist=dto.edgelist;
  vector <double> trianglelist=dto.trianglelist;
  vector <double> pointlist=dto.pointlist;
  vector <double> pointattributelist=dto.pointattributelist;
  vector <int> pointmarkerlist=dto.pointmarkerlist;
  vector <double> regionlist=dto.regionlist;
  vector <int> segmentlist=dto.segmentlist;
  vector <int> segmentmarkerlist=dto.segmentmarkerlist;
  double maximumTriangleArea=dto.maximumTriangleArea;
  Tricall_Utility tri;
  struct triangulateio in, mid, vorout;
  /* Define input points. */

  in.numberofpoints = dto.numberofpoints;;
  in.numberofpointattributes = dto.numberofpointattributes;;
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
  
  in.numberofsegments = 0;
  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  for(int iTer=0;iTer < (in.numberofsegments * 2 );iTer++){
  in.segmentlist[iTer]=segmentlist[iTer];
  }
  
  in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(REAL));
  for(int iTer=0;iTer < (in.numberofsegments);iTer++){
  in.segmentmarkerlist[iTer]=segmentmarkerlist[iTer];
  }
  
  in.numberofholes = dto.numberofholes;
  
  in.numberofregions = dto.numberofregions;
  in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
  for(int iTer=0;iTer < (in.numberofregions * 4 );iTer++){
  in.regionlist[iTer]=regionlist[iTer];
  }          /* Area constraint that will not be used. */
  
  printf("IIIIIIIIIInput point set:\n\n");
  //report(&in, 1, 0, 0, 0, 0, 0);

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
    
  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), number everything from  */
  /*   zero (z), assign a regional attribute to each element (A), and  */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  std::ostringstream oss; // string as stream
  //oss << "a" << maximumTriangleArea << "AenYPQ"; // write to string stream
  oss << "a" << maximumTriangleArea << options_triangle_second_call; // write to string stream
  string flagsStr = oss.str(); // get string out of stream
  char * flagsStrChars = new char [flagsStr.length()];
  std::strcpy (flagsStrChars, flagsStr.c_str());
  char *filename=dto.filenamestr;
  triangulate(flagsStrChars, &in, &mid, (struct triangulateio *) NULL, filename, 1);
  
  //Fill data to be sent to tricall pre mesh gen - Start
  //Fill points details
  pointlist.erase(pointlist.begin(),pointlist.end());
  for (int i = 0; i < mid.numberofpoints; i++)
  {
  pointlist.push_back(mid.pointlist[i * 2 + 0]);
  pointlist.push_back(mid.pointlist[i * 2 + 1]);
  }
  dto.numberofpoints=mid.numberofpoints;
  dto.pointlist=pointlist;
  
  //Fill edge lists
  //Wont work with 'z' flag
  for (int i = 0; i < mid.numberofedges; i++)
  {
  edgelist.push_back((mid.edgelist[i * 2 + 0])-1);
  edgelist.push_back((mid.edgelist[i * 2 + 1])-1);
  }
  dto.numberofedges=mid.numberofedges;
  dto.edgelist=edgelist;
  
  //Fill element lists
  //Wont work with 'z' flag
  for (int i = 0; i < mid.numberoftriangles; i++)
  {
  for (int j = 0; j < mid.numberofcorners; j++)
  {
    trianglelist.push_back((mid.trianglelist[i * mid.numberofcorners + j])-1);
  }
  }
  dto.numberoftriangles=mid.numberoftriangles;
  dto.trianglelist=trianglelist;
  dto.numberofcorners=mid.numberofcorners;
  
  tri.exportnodes(&mid,1,filename);
  tri.exportedges(&mid,1,filename);
  tri.exportelements(&mid,filename);
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
  free(mid.neighborlist);
  free(mid.segmentlist);
  free(mid.segmentmarkerlist);
  free(mid.edgelist);
  free(mid.edgemarkerlist);
  std::cout << "Done From generateMeshWithFixedBoundaries for " << dto.filenamestr << "\n" ;

}

//This is currently used
void Tricall_Utility::generateTriangleMesh(Data_Transfer_Object &dto)
{
  std::cout << "Entered Tricall_Utility::generateTriangleMesh for " << dto.filenamestr <<"\n" ;
  vector <double> edgelist=dto.edgelist;
  vector <double> trianglelist=dto.trianglelist;
  vector <double> pointlist=dto.pointlist;
  vector <double> pointattributelist=dto.pointattributelist;
  vector <int> pointmarkerlist=dto.pointmarkerlist;
  vector <double> regionlist=dto.regionlist;
  vector <int> segmentlist=dto.segmentlist;
  vector <int> segmentmarkerlist=dto.segmentmarkerlist;
  double maximumTriangleArea=dto.maximumTriangleArea;
  
  /**printf("Input point set Before triangulate:\n\n");
  for(int i=0; i< edgelist.size();i++)
  {
  cout << edgelist[i] << "~~";
  }
  cout << "\n";
  for(int i=0; i< trianglelist.size();i++)
  {
  cout << trianglelist[i] << "!!";
  }
  cout << "\n";
  for(int i=0; i< pointlist.size();i++)
  {
  cout << pointlist[i] << "@@";
  }
  cout << "\n";
  for(int i=0; i< pointattributelist.size();i++)
  {
  cout << pointattributelist[i] << "##";
  }
  cout << "\n";
  for(int i=0; i< pointmarkerlist.size();i++)
  {
  cout << pointmarkerlist[i] << "$$";
  }
  cout << "\n";
  for(int i=0; i< regionlist.size();i++)
  {
  cout << regionlist[i] << "%%";
  }
  cout << "\n";
  for(int i=0; i< segmentlist.size();i++)
  {
  cout << segmentlist[i] << "^^";
  }
  cout << "\n";
  for(int i=0; i< segmentmarkerlist.size();i++)
  {
  cout << segmentmarkerlist[i] << "&&";
  }
  cout << "\n";**/
  
  
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

  /**printf("Input point set Before triangulate:\n\n");
  cout << pointlist.size() << "\n";
  for(int i=0; i< pointlist.size();i++)
  {
  cout << pointlist[i] << "::";
  }
  for(int i=0; i< segmentlist.size();i++)
  {
  cout << segmentlist[i] << "##";
  }
  cout << "\n";**/
  
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

  triangulate(options_triangle_first_call, &in, &mid, &vorout, "test", 0);
 
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
  //triangulate("enpraDPQ", &mid, &out, (struct triangulateio *) NULL, filename, 1);
  char *meshParam=dto.meshCommandLine;
  cout << meshParam << "" << filename << "\n";
  triangulate(meshParam, &mid, &out, (struct triangulateio *) NULL, "test", 1);
  
  //printf("Refined triangulation:\n\n");
  //tri.report(&out, 1, 1, 1, 0, 1, 0);
  
  //Fill data to be sent to tricall pre mesh gen - Start
  //Fill points details
  pointlist.erase(pointlist.begin(),pointlist.end());
  for (int i = 0; i < out.numberofpoints; i++)
  {
  pointlist.push_back(out.pointlist[i * 2 + 0]);
  pointlist.push_back(out.pointlist[i * 2 + 1]);
  }
  dto.numberofpoints=out.numberofpoints;
  dto.pointlist=pointlist;
  
  //TODO File point attribute and marker lists
  
  //Fill edge lists
  //Wont work with 'z' flag
  for (int i = 0; i < out.numberofedges; i++)
  {
  edgelist.push_back((out.edgelist[i * 2 + 0])-1);
  edgelist.push_back((out.edgelist[i * 2 + 1])-1);
  }
  dto.numberofedges=out.numberofedges;
  dto.edgelist=edgelist;
  
  //Fill element lists
  //Wont work with 'z' flag
  for (int i = 0; i < out.numberoftriangles; i++)
  {
  for (int j = 0; j < out.numberofcorners; j++)
  {
    trianglelist.push_back((out.trianglelist[i * out.numberofcorners + j])-1);
  }
  }
  dto.numberoftriangles=out.numberoftriangles;
  dto.trianglelist=trianglelist;
  dto.numberofcorners=out.numberofcorners;
  
  //Fill data to be sent to tricall pre mesh gen - End
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
  
  std::cout << "Exited Tricall_Utility::generateTriangleMesh for " << dto.filenamestr << "\n" ;
}

//Wont work if 'z' flag is used
void Tricall_Utility::exportelements(struct triangulateio *io,char *filename)
{
  //printf("Entered export elements\n");
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
    
    //printf("Triangle %4d points:", i);
    fprintf(outfile, "%4ld    ", i);
    for (int j = 0; j < io->numberofcorners; j++) {
      //printf("  %4d", (io->trianglelist[i * io->numberofcorners + j])-1);
      fprintf(outfile, "%4d  ", (io->trianglelist[i * io->numberofcorners + j])-1);
    }
    if (io->numberoftriangleattributes > 0) {
      //printf("   attributes");
    }
    for (int j = 0; j < io->numberoftriangleattributes; j++) {
      //printf("  %.6g", io->triangleattributelist[i *
                     //io->numberoftriangleattributes + j]);
      fprintf(outfile, "%.17g  ", io->triangleattributelist[i *
                     io->numberoftriangleattributes + j]);								 
    }
    //printf("\n");
    fprintf(outfile, "\n");		
  }
  //printf("\n");
  fclose(outfile);
}

//Wont work if 'z' flag is used
void Tricall_Utility::exportedges(struct triangulateio *io, int markers, char *filename)
{
  //printf("Entered export edges\n");
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
  //printf("Mid writing edge files\n");
  
  for (int i = 0; i < io->numberofedges; i++) {
      //printf("Edge %4d points:", i);
    fprintf(outfile, "%4ld   ", i);
      for (int j = 0; j < 2; j++) {
        //printf("  %4d", (io->edgelist[i * 2 + j])-1);
    fprintf(outfile, "%d  ", (io->edgelist[i * 2 + j])-1);
      }
      
      if (markers) {
        //printf("   marker %d\n", io->edgemarkerlist[i]);
    fprintf(outfile, "%d\n", io->edgemarkerlist[i]);
      } else {
        //printf("\n");
    fprintf(outfile, "\n");
      }
    }
    //printf("\n");
  printf("\n");
  fclose(outfile);
  //printf("Exiting exportnodes\n");
}

void Tricall_Utility::exportnodes(struct triangulateio *io, int markers, char *filename)
{
  int i, j;

  //printf("Entered export nodes\n");
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
    //printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      //printf("  %.6g", io->pointlist[i * 2 + j]);
    }
  
  /* Vertex number, x and y coordinates. */
    fprintf(outfile, "%4d    %.17g  %.17g", i, io->pointlist[i * 2 + 0],
              io->pointlist[i * 2 + 1]);
  
    if (io->numberofpointattributes > 0) {
      //printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      //printf("  %.6g",
             //io->pointattributelist[i * io->numberofpointattributes + j]);
    /* Write an attribute. */
      fprintf(outfile, "  %.17g", io->pointattributelist[i * io->numberofpointattributes + j]);		 
    }
    if (markers) {
      //printf("   marker %d\n", io->pointmarkerlist[i]);
    /* Write the boundary marker. */
      fprintf(outfile, "    %d\n", io->pointmarkerlist[i]);
    } else {
      //printf("\n");
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

void Tricall_Utility::exportnodesfrmArray(vector<double>& pointlist,char *filename,int numberofpoints)
{
        
    FILE *outfile;
    
    char *extension="_st.node";
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
    
    fprintf(outfile, "%d  %d  %d  %d\n", numberofpoints, 2,
        0, 0);
    
    for(int iPoint=0;iPoint < numberofpoints ;iPoint++)
    {
    double coord1=pointlist[(iPoint*2)];
    double coord2=pointlist[(iPoint*2)+1];
    fprintf(outfile, "%4d    %.17g  %.17g", iPoint, coord1,coord2);
    fprintf(outfile, "\n");
    }
  
    fclose(outfile);
    
    
    
    
    
    
}	  







void Tricall_Utility::exporttrianglesfrmArray(int *triangle_ptr,char *filename,int numberoftriangles,int numberofcorners)
{
  //printf("Entered export elements\n");
  FILE *outfile;

  //Start concatenation
  //Concatenation two char* strings
  //To dynamically create the filename with ".node" appended 
  //to the end of the filename
  char *extension="_st.ele";
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
  fprintf(outfile, "%d  %d  %d\n", numberoftriangles,
      numberofcorners, 0);
  
  for (int i = 0; i < numberoftriangles; i++) {
    
    if((*(triangle_ptr + (i*(numberofcorners+1)) + 3)) == 1)
    
    {
    
      fprintf(outfile, "%4ld    ", i);
      for (int j = 0; j < numberofcorners; j++) {
        fprintf(outfile, "%4d  ", (*(triangle_ptr + (i*(numberofcorners+1)) + j)));
      }
      
      fprintf(outfile, "\n");	
    }
    
  }
  
  fclose(outfile);
}

void Tricall_Utility::exportedgesfrmArray(int *edgearr,int *newedgearr,char *filename,int numberofedges,int numberofedgesnew)
{
  //printf("Entered export edges\n");
  FILE *outfile;

  //Start concatenation
  //Concatenation two char* strings
  //To dynamically create the filename with ".node" appended 
  //to the end of the filename
  char *extension="_st.edge";
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
  
  int totalnumberofedges=numberofedges+numberofedgesnew;
  
  /* Number of edges, number of boundary markers (zero or one). */
  fprintf(outfile, "%d  %d\n", totalnumberofedges, 0);
  
  
  for (int iEdgeIter=0;iEdgeIter < numberofedges;iEdgeIter++) {
      fprintf(outfile, "%4ld   ", iEdgeIter);
    int node1=(*(edgearr + (iEdgeIter*2)));
    int node2=(*(edgearr + (iEdgeIter*2) + 1));
      fprintf(outfile, "%d  %d  ", node1,node2);
      fprintf(outfile, "\n");      
    }
  
  for (int iEdgeIter=0;iEdgeIter < numberofedgesnew;iEdgeIter++) {
      fprintf(outfile, "%4ld   ", numberofedges+iEdgeIter);
    int node1=(*(newedgearr + (iEdgeIter*2)));
    int node2=(*(newedgearr + (iEdgeIter*2) + 1));
      fprintf(outfile, "%d  %d  ", node1,node2);
      fprintf(outfile, "\n");      
    }
  
  fclose(outfile);
  
}




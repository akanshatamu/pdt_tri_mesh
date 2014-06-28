#include <vector>
#include <string>
#include <cmath>
#include <map>


using std::vector;
using std::string;
using std::map;

class Data_Transfer_Object
{
public:
	Data_Transfer_Object(void);
	~Data_Transfer_Object(void);
	
	int numberofpoints;
	int numberofpointattributes;
	vector<double> pointlist;
	vector<double> pointattributelist;
	vector<int> pointmarkerlist;
	int numberofsegments;
	int numberofholes;
	int numberofregions;
	vector<double> regionlist;
	double maximumTriangleArea;
	char *filenamestr;
	vector<int> segmentlist;
	vector<int> segmentmarkerlist;
	
	
	
};


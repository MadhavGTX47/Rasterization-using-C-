#ifndef stuff_H
#define stuff_H
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>


using namespace std;
using namespace Eigen;

class View{			//Has all the basic values
	public:
		double b[3];
		double f[8];
		Vector3d from;
		Vector3d up;
		Vector3d at;
		double angle;
		double hither;
		int resx,resy;
		

};

class Triangleraster  {   //Extra structure to store transfomred points and colors
    
    public:
  
        Vector3d points[3]; 
        Vector3d norms[3]; 
        Vector3d color[3]; 

};

struct fragment { 		//Structure for storing fragments
    double zval;     
    Vector3d color;     
    struct fragment *next;
};

int rasterizer();

vector<Matrix<double,9, 3>> tri;
vector<Triangleraster> trimage;
vector<Vector3d> light;
View view;  //View Object which stores our vectors and other important data

#endif

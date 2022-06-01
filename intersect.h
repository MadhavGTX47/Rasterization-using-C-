#ifndef intersect_H
#define intersect_H
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include "stuff.h"

using namespace std;
using namespace Eigen;

Vector3d shade(Matrix<double,9, 3> tri,int row)
{   
        Vector3d a = tri.row(0);
	    Vector3d b = tri.row(2);
	    Vector3d c = tri.row(4);
       
        Vector3d ret = {0,0,0};
        Vector3d vertex = tri.row(row);           
        Vector3d normal;


        double intensity = 1/sqrt(light.size());

	    Vector3d view_direction =view.from-vertex;     //Calculating View direction
		view_direction= view_direction.normalized();	

        Vector3d vector_avg ;	

        for(int j = 0; j < light.size(); j++){
    			
		Vector3d light_direction = light[j]-vertex;    //calculating light direction
		light_direction= light_direction.normalized();
        
        if(tri(8,2)==1)
        {
         normal =((b-a).cross(c-a)).normalized(); 

        }
        else{
         normal=(tri.row(row+1));
        }

        vector_avg = (view_direction + light_direction)/(view_direction + light_direction).norm();   //calculating view direction
        
			

        double diffuse = max(0.0, normal.dot(light_direction));  //Shading colors
        double specular = max(0.0,(normal.dot(vector_avg)));
		specular = tri(7,1)*pow(specular,tri(7,2));

     

            ret[0]+=(tri(7,0) * tri(6,0) * diffuse) * intensity;  
            ret[1]+=(tri(7,0) * tri(6,1) * diffuse) * intensity;
            ret[2]+=(tri(7,0) * tri(6,2) * diffuse) * intensity;  

            
            ret[0]+=specular * intensity;
            ret[1]+=specular * intensity;
            ret[2]+=specular * intensity;
        }  

      if(ret[0]>1) 									//Limit the colors to be in (0,1)
		{ret[0]=1; }
		if(ret[1]>1)
		{ret[1]=1; }
		if(ret[2]>1)
		{ret[2]=1; } 
     
      if(ret[0]<0) 									//Limit the colors to be in (0,1)
		{ret[0]=0; }
		if(ret[1]<0)
		{ret[1]=0; }
		if(ret[2]<0)
		{ret[2]=0; }    
        
        return ret; 
}


fragment* fragz(int x, int y,Matrix<double,9, 3>  &t, Triangleraster &tr) {
    
    fragment *frag = NULL;                      //Interpolating colors

    double x0 = tr.points[0][0];
    double y0 = tr.points[0][1];
    double x1 = tr.points[1][0];
    double y1 = tr.points[1][1];
    double x2 = tr.points[2][0]; 
    double y2 = tr.points[2][1];
   
    double beta = ((y2-y0)*x + (x0-x2)*y + x2*y0 - x0*y2) / ((y2-y0)*x1 + (x0-x2)*y1 + x2*y0 - x0*y2);
    double gamma = ((y0-y1)*x + (x1-x0)*y + x0*y1 - x1*y0) / ((y0-y1)*x2 + (x1-x0)*y2 + x0*y1 - x1*y0);
    double alpha=1-beta-gamma;
   
    
    if (alpha>0 && beta>0 && gamma>0) {             //Assigining colors for the fragments along with zvalue and transperancy
        frag = new fragment; 
        frag->color = alpha*tr.color[0]+ beta*tr.color[1] + gamma*tr.color[2];
        frag->zval = alpha*tr.points[0][2]+ beta*tr.points[1][2]+ gamma*tr.points[2][2];
        frag->next = NULL;
    }

    return frag;
}
/*
fragment* fragl(int x, int y,Matrix<double,9, 3>  &tri, Triangleraster &tr) {
    
    fragment *frag = NULL;                      //Interpolating colors
        Vector3d a = {tr.points[0][0],tr.points[0][1],tr.points[0][2]};
	    Vector3d b = {tr.points[1][0],tr.points[1][1],tr.points[1][2]};
	    Vector3d c = {tr.points[2][0],tr.points[2][1],tr.points[2][2]};
       
        Vector3d ret = {0,0,0};
for(int i;i<3;i++)
{
        Vector3d vertex ={tr.points[i][0],tr.points[i][1],tr.points[i][2]};          
        Vector3d normal;


        double intensity = 1/sqrt(light.size());

	    Vector3d view_direction =view.from-vertex;     //Calculating View direction
		view_direction= view_direction.normalized();	

        Vector3d vector_avg ;	

        for(int j = 0; j < light.size(); j++){
    			
		Vector3d light_direction = light[j]-vertex;    //calculating light direction
		light_direction= light_direction.normalized();
        
        if(tri(8,2)==1)
        {
         normal =((b-a).cross(c-a)).normalized(); 

        }
        else{
         normal=(tri.row(i*2+1));
        }

        vector_avg = (view_direction + light_direction)/(view_direction + light_direction).norm();   //calculating view direction
        
			

        double diffuse = max(0.0, normal.dot(light_direction));  //Shading colors
        double specular = max(0.0,(normal.dot(vector_avg)));
		specular = tri(7,1)*pow(specular,tri(7,2));

     

            ret[0]+=(tri(7,0) * tri(6,0) * diffuse) * intensity;  
            ret[1]+=(tri(7,0) * tri(6,1) * diffuse) * intensity;
            ret[2]+=(tri(7,0) * tri(6,2) * diffuse) * intensity;  

            
            ret[0]+=specular * intensity;
            ret[1]+=specular * intensity;
            ret[2]+=specular * intensity;
        }  

      if(ret[0]>1) 									//Limit the colors to be in (0,1)
		{ret[0]=1; }
		if(ret[1]>1)
		{ret[1]=1; }
		if(ret[2]>1)
		{ret[2]=1; } 
     
      if(ret[0]<0) 									//Limit the colors to be in (0,1)
		{ret[0]=0; }
		if(ret[1]<0)
		{ret[1]=0; }
		if(ret[2]<0)
		{ret[2]=0; } 

    tr.color[i]=ret;
}
    double x0 = tr.points[0][0];
    double y0 = tr.points[0][1];
    double x1 = tr.points[1][0];
    double y1 = tr.points[1][1];
    double x2 = tr.points[2][0]; 
    double y2 = tr.points[2][1];
   
    double beta = ((y2-y0)*x + (x0-x2)*y + x2*y0 - x0*y2) / ((y2-y0)*x1 + (x0-x2)*y1 + x2*y0 - x0*y2);
    double gamma = ((y0-y1)*x + (x1-x0)*y + x0*y1 - x1*y0) / ((y0-y1)*x2 + (x1-x0)*y2 + x0*y1 - x1*y0);
    double alpha=1-beta-gamma;
   
    
    if (alpha>0 && beta>0 && gamma>0) {             //Assigining colors for the fragments along with zvalue and transperancy
        frag = new fragment; 
        frag->color = alpha*tr.color[0]+ beta*tr.color[1] + gamma*tr.color[2];
        frag->zval = alpha*tr.points[0][2]+ beta*tr.points[1][2]+ gamma*tr.points[2][2];
        frag->next = NULL;
    }

    return frag;
}
*/
#endif

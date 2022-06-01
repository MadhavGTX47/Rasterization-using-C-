#include <string>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>
#include <math.h>       
#include <Eigen/Geometry>
#include "stuff.h"
#include "intersect.h"
using namespace Eigen;
using namespace std;



double trans = -1;									 

int main(int argc, char* argv[])
{	 

 int c;
  while ((c = getopt (argc, argv, "sft:")) != -1)  //Handling command line flags
    switch (c)
      {
      case 's':
       cout<<"s";
        break;
      case 'f':
        cout<<"f";
        break;	
      case 't':
        trans = atof(optarg);;
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }


	
	if(argv[1])
	  cout<<"Generating the output.ppm";		//To take the Input file name
	else
	  cout << "Give proper arguments" << endl;
	

    Matrix<double,9, 3> triangle; 		  		 // Temporary variable to store Traingle data while entering into vector
	RowVector3d l;								 // Temporary variable to store light data while entering into vector			
    fstream nff(argv[optind++]);           		 // Making an fstream for the nff file 
	string line;               			  		 // to get each line from nff file
	vector<string> tokens;     			  		 // vector to store the tokens
	string temp;                          		 // temporay string used in tokenization
	
	
	
	while (std::getline(nff, line)){        	 //to read the nff file line by line
	
		istringstream iss(line);
		stringstream li(line);          		 // storing the line in li of stringstream
	
	while(getline(li, temp, ' ')){       		 // tokenization with white spaces as delimter
        tokens.push_back(temp);           		 // storing the tokens in the token vector
    }

	}

	
	for(int i=0; i < tokens.size(); i++){  		//Iterating through tokens
	   
	   
		if(tokens[i]=="b")					    //Storing the background color
		{
			view.b[0]=stod(tokens[++i]);
			view.b[1]=stod(tokens[++i]);
			view.b[2]=stod(tokens[++i]);
		}
		
		if(tokens[i]=="v")					    //Storing the view
		{
			i++;			
			if(tokens[i]=="from")				//Storing the from coordinate
			{
			view.from<<stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]); 
			}
			i++;
			if(tokens[i]=="at")					//Storing the At coordinate
			{
			view.at<<stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]);  
			}
			i++;
			if(tokens[i]=="up")					//Storing the Up coordinate
			{
			view.up<<stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]);  
			}
			i++;
			if(tokens[i]=="angle")				//Storing the Angle
			{
			view.angle=stod(tokens[++i]);
			}
			i++;
			if(tokens[i]=="hither")				//Storing the Hither
			{
			view.hither=stod(tokens[++i]);  
			}
			i++;
			if(tokens[i]=="resolution")			//Storing the Resolution
			{
			view.resx=stoi(tokens[++i]); view.resy=stoi(tokens[++i]);
			}
		
		}
		
		if(tokens[i]=="l")						//Storing the lights
			{
			l<< stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]);
			light.push_back(l);
			}	



		if(tokens[i]=="f")					    //Storing the fill color
		{
			view.f[0]=stod(tokens[++i]);
			view.f[1]=stod(tokens[++i]);
			view.f[2]=stod(tokens[++i]);
			view.f[3]=stod(tokens[++i]);
			view.f[4]=stod(tokens[++i]);
			view.f[5]=stod(tokens[++i]);
			view.f[6]=stod(tokens[++i]);
			view.f[7]=stod(tokens[++i]);

		}
		
		
		if(tokens[i]=="p")					    //Storing the polygons
		{	
			double ver;							//temp vector
			int v = stoi(tokens[++i]) ;         //Converting token which is in string to integer for comparison
			vector<float> tr;

			if(v>3)								//To handle polygons with more than 3 vertices
			{ 
				for(int tx=0;tx<v*3;tx++)		//Storing all the vertices in temporary vector to split and store
				{ 
				ver = stod(tokens[++i]);
				tr.push_back(ver);
				}
								
				v=v-2;
				int t=0;
				for(double co=0;co < v;co++)
				{
				triangle << tr[0],tr[1],tr[2],0,0,0,tr[t+3],tr[t+4],tr[t+5],0,0,0,tr[t+6],tr[t+7],tr[t+8],0,0,0,view.f[0],view.f[1],view.f[2],view.f[3],view.f[4],view.f[5],view.f[6],view.f[7],1;
				tri.push_back(triangle);
				t=t+3;							//Here we saved our fill color in last row of the matrix
			    }
			}
			
			else
			{   								// To handle polygons with 3 vertices
				triangle << stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,0,stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,0,stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),0,0,0,view.f[0],view.f[1],view.f[2],view.f[3],view.f[4],view.f[5],view.f[6],view.f[7],0; 
				tri.push_back(triangle);		//Here we saved our fill color in last row of the matrix
				
			}
			
		}
	
		if(tokens[i]=="pp")					    //Storing the polygon patches
		{	
			double ver;							//temp vector
			int v = stoi(tokens[++i]) ;         //Converting token which is in string to integer for comparison
			vector<float> tr;

			if(v>3)								//To handle polygon patch with more than 3 vertices
			{ 
				cout<<"Currently not being handled";		//Storing all the vertices in temporary vector to split and store
				
			}
			
			else
			{   								// To handle polygons patches with 3 vertices
				triangle << stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),stod(tokens[++i]),view.f[0],view.f[1],view.f[2],view.f[3],view.f[4],view.f[5],view.f[6],view.f[7],0; 
				tri.push_back(triangle);		//Here we saved our fill color in last row of the matrix
				
			}
			
		}	
		
	}

	rasterizer();      						   //Calling our raytracer as we have finished parsing

return 0;

}


int rasterizer() {		
	
	fragment ***frags = new fragment**[view.resy];		  //Initializing fragments
    	for (int j=0; j<view.resy; j++)
		{      
        
		frags[j] = new fragment*[view.resx];
       
	    for (int i=0; i<view.resx; i++) 
		{   
          
		    frags[j][i] = NULL;							//Setting all fragments to null initially
        }
   		}



	unsigned char pixels[view.resx][view.resy][3];	      // pixels which is to store output

	view.up=-1*view.up;
	Vector3d w = (view.from-view.at).normalized(); 		  // Calculating u,v,w
	Vector3d u = -((view.up.cross(w)).normalized());
	Vector3d v = -((w.cross(u)).normalized());
	
	int nx = view.resx;
    int ny = view.resy;
	double AspectRatio = view.resx / view.resy;

	double h = (tan ((view.angle * M_PI / 180.0) /2));	  // Assigning and calcualting h,r,l,t,b
    double right =h;
    double left = -1.0 * right;
    double top = right / AspectRatio;
    double bottom = left / AspectRatio;

	double near =-1* view.hither;
    double far = 1* 1000 * near;

    Matrix4d Mvp, Mper;									//Calucalting Mvp,Mper,Mcam
    MatrixXd Mcam; 
	Mcam.resize(3,4);

    Mvp << (nx/2.0), 0, 0, ((nx-1)/2.0),0, ((ny/2.0)), 0, ((ny-1)/2.0),0, 0, 1, 0,0, 0, 0, 1;

    
    Mper << (2.0*near)/(right-left), 0, (left+right)/(left-right), 0,0, (2.0*near)/(top-bottom), (bottom+top)/(bottom-top), 0,0, 0, (far+near)/(near-far), (2.0*far*near)/(far-near),0, 0, 1, 0;
	
	
    Mcam << u, v, w, view.from;
    Mcam.conservativeResize(4, 4);
    Mcam.row(3) = Vector4d(0,0,0,1);
    Mcam = Mcam.inverse();

	Matrix4d M = Mvp * Mper * Mcam;                    	 //Calculating M


	vector<Triangleraster> trimage(tri.size());

	for (long int i=0; i<tri.size(); i++) {
        Matrix<double,9, 3> tria = tri[i];  

        for (int j=0; j<3; j++) {

           	trimage[i].color[j] = shade(tria, j*2);     //Calculating colors at vertex

            
            Vector4d temp;
			Vector3d p=tria.row(j*2);    				// transforming vertices from world space to image space
       		temp << p, 1.0;
            temp = M * temp;
            temp = temp / temp[3];  					//  homogenous divide 
            trimage[i].points[j] = {temp[0], temp[1], temp[2]};
			
			
        }
    }



	
	for (unsigned int i=0; i<tri.size(); i++) {
      
	    Matrix<double,9, 3> &trif = tri[i];  
		Triangleraster	&t=trimage[i];
 
        double minx=t.points[0][0];					//Calculating bounding box
		double maxx=minx;
		double miny=t.points[0][1];
		double maxy=miny;

        for (int n=1; n<3; n++) {	

            if (t.points[n][0] < minx)
			{
                minx = t.points[n][0];
			}	
            if (t.points[n][0] > maxx)
			{
                maxx = t.points[n][0];
			}	
            if (t.points[n][1] < miny)
			{
                miny = t.points[n][1];
			}	
            if (t.points[n][1] > maxy)
			{
                maxy = t.points[n][1];
			}	
        }

        minx = floor(minx);
        miny = floor(miny);
        maxx = ceil(maxx);
        maxy = ceil(maxy);


		

        for (int y=(int) miny; y<=(int) maxy; y++) {		//Checking bounding boxes
            for (int x= (int)minx; x<= (int)maxx; x++) {
               
                if (!(0<=x && x<view.resx && 0<=y && y<view.resy)) {
                    continue;
                }

				fragment *fragt = fragz(x, y, trif,t);    //Storing fragments in a linked list
				if (fragt == NULL) 
                    continue;
               
                if (frags[y][x] == NULL) 
				{
                    frags[y][x] = fragt;
                }
                else 
				{
                    fragt->next = frags[y][x]; 			  //Attaching fragments in as a list
                    frags[y][x] = fragt;
                }
				
            }
        }
    
	}
if(trans==-1){
	for (int iy=0; iy<view.resy; iy++) {      
        for (int ix=0; ix<view.resx; ix++) {   //Calcualting fragments
            
            Vector3d color ={view.b[0],view.b[1],view.b[2]};

            if (frags[iy][ix] != NULL) {
            
                fragment *best=frags[iy][ix];
				fragment *te=best;
                while (te != NULL) {
                    
                    if (te->zval > best->zval) {	//Selecting the nearest Z value
                        best = te;
                    }
                    te = te->next;
                }
                color = best->color;
            }
                color[0] = min(max(color[0], 0.0), 1.0); 
				color[1] = min(max(color[1], 0.0), 1.0);
				color[2] = min(max(color[2], 0.0), 1.0);
				
				pixels[iy][ix][0]= color[0]*255;        //Assiging colors to the pixels
				pixels[iy][ix][1]= color[1]*255;
				pixels[iy][ix][2]= color[2]*255;
				
            }
        }
}
else{
	for (int iy=0; iy<view.resy; iy++) {          //for transperent teapot
        for (int ix=0; ix<view.resx; ix++) {  
            
            Vector3d color ={view.b[0],view.b[1],view.b[2]};
			
            if (frags[iy][ix] != NULL) { 
            
                fragment *te=frags[iy][ix];
				fragment *fs=te;
				fragment *s=te;

			while (fs != NULL) {				//Sorting by  Z-value the fragments

			fragment *best=fs;
			s=fs;
			Vector3d fcolor;
			while (s != NULL) {

				if (s->zval < best->zval) {			//basic selection sort
					double tempz;
					Vector3d tempc;
					tempz=best->zval;
					tempc=best->color;
					best->zval=s->zval;
					best->color=s->color;
					s->color=tempc;
					s->zval=tempz;
				}
                    s = s->next;
                }

				te=best;
				
				fcolor = te->color;                //Transperancy color addition
                color[0]=(trans*fcolor[0])+(1.0-trans)*color[0];
				color[1]=(trans*fcolor[1])+(1.0-trans)*color[1];
				color[2]=(trans*fcolor[2])+(1.0-trans)*color[2];
				fs=fs->next;
				}
				
				
            }
     
				color[0] = min(max(color[0], 0.0), 1.0);  // Storing colors in the pixels array
				color[1] = min(max(color[1], 0.0), 1.0);
				color[2] = min(max(color[2], 0.0), 1.0);
		
				pixels[iy][ix][0]= color[0]*255;
				pixels[iy][ix][1]= color[1]*255;
				pixels[iy][ix][2]= color[2]*255;
			
				
            }
			
        }
}



	
		cout<<endl<<endl<<"Writing...";
		FILE *f = fopen("output.ppm","wb");		   //Saving the output to a ppm file
		fprintf(f, "P6\n%d %d\n%d\n", view.resx, view.resy, 255);
		fwrite(pixels, 1, view.resx*view.resy*3, f);
		fclose(f);	
		cout<<endl<<endl<<"done. ";



return 0;
}		





		




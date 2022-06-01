raster : raster.o 
	g++ raster.o -o raster

raster.o : raster.cpp 
	g++ -I/usr/include/eigen3/ -O3 -c raster.cpp


clean:
	rm *.o *.ppm *.exe output out
	 

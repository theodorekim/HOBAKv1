SHELL := /bin/bash -e

mac_omp : 
	echo -e "\n==== Building for OSX with OpenMP ====\n";cd projects; cp -f include_top.mac_omp include_top.mk; cd simulateScene; make depend; make clean; make -j4; cd ../../;cd bin;./simulateScene

mac : 
	echo -e "\n==== Building for OSX ====\n";cd projects; cp -f include_top.mac include_top.mk; cd simulateScene; make depend; make clean; make -j4; cd ../../;cd bin;./simulateScene

linux : 
	echo -e "\n==== Building for Linux ====\n";cd projects; cp -f include_top.linux include_top.mk; cd simulateScene; make depend; make clean; make -j4; cd ../../;cd bin;./simulateScene

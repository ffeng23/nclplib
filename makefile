install: clplib

objects:
	g++ -std=c++11 -Wall -O3 -fPIC -c *.cpp

clplib: objects
	g++ -shared -Wall -std=c++11 -O3 -Wl,-soname,libCLPLib.so.1 -o libCLPLib.so.1.0 *.o

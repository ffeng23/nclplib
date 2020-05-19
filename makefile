LIBCLP= libCLPLib.so.1.0

all: $(LIBCLP) 

.PHONY: clean all install objects clplib

clean:
	rm -fr *.o *~ core
	#cd $(ACCESSDIR) && rm -fr *.o *~ core; #ls ;   #pwd ; # note makefile is a script and each line (of commands) is doing its own sub-precess (&& is trying to stop if the "cd" command fails.


install: $(LIBCLP)

objects: 
	g++ -std=c++11 -Wall -O3 -fPIC -c *.cpp

clplib: objects
	g++ -shared -Wall -std=c++11 -O3 -Wl,-soname,libCLPLib.so.1 -o libCLPLib.so.1.0 *.o

$(LIBCLP): clplib
	
#Make file for doing compiling of nclplib c++ library
#this depends on the conda environment "source activate FengPython3.7" on theano (-L/home/ffeng/miniconda3/lib/)
#	
# but use the system settings on P52S.
#need to make sure the compiler is correct based on environment
#for conda evironment (doing nothing, since codan will set up the correct g++ compiler using GXX
#
#for using systemwise g++ <----------
GXX=g++

#this is try to use conda env variblbe. <------ uncomment on P52s
#LD=${LDFLAGS} -L/home/ffeng/miniconda3/lib/
LD=
LIBCLP= libCLPLib.so.1.0

all: $(LIBCLP) 

.PHONY: clean all install objects clplib

clean:
	rm -fr *.o *~ core
	#cd $(ACCESSDIR) && rm -fr *.o *~ core; #ls ;   #pwd ; # note makefile is a script and each line (of commands) is doing its own sub-precess (&& is trying to stop if the "cd" command fails.


install: $(LIBCLP)

objects: 
	${GXX} -std=c++11 -Wall -O3 -fPIC -m64 -c *.cpp

clplib: objects
	${GXX} -shared -Wall -std=c++11 -O3 -m64 ${LD} -Wl,-soname,libCLPLib.so.1  -o libCLPLib.so.1.0 *.o

$(LIBCLP): clplib


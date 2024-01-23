#include <iostream>
#include <dlfcn.h>

typedef int (*moo)();
extern "C"
{
void print_line(const char *str);
}



int main()
{
  std::cout << "Hello World!" << std::endl;

  //print_line("this is the line print from lib");


  void* myso = dlopen("/usr/lib/libCLPLib.so", RTLD_NOW );

  if (!myso){
    std::cout << "Failed to load lib " << dlerror() << std::endl;
    return 1;
  }
  else 
  {

	//print_line("this is the line print from lib");
	std::cout << "open lib ok\n";

   }


  dlclose(myso);

  return 0;
}

#include <chrono>
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>

#define START(x) auto x##start = std::chrono::steady_clock::now()                                         
#define END(x)   auto x##end   = std::chrono::steady_clock::now()                                         
#define TIME(x)  auto x##time  = std::chrono::duration_cast<std::chrono::milliseconds>(x##end-x##start).count()
//#define PRINT(x) std::cout << #x << "(), " << x##time << std::endl
//#define PRINT(x) std::string x##print(#x + "(), " + x##time); std::cout << x##print << std::endl

#define PRINT(x) std::stringstream x##_ss; x##_ss << #x << "(), " << x##time << std::endl;\
                 std::string x##_str = x##_ss.str(); \
                 std::cout << x##_str;

#define COMP(X)  END(X); TIME(X); PRINT(X);

#define CSTART(x) \
    struct timespec x##start;\
    clock_gettime(CLOCK_REALTIME, &x##start); \

#define CEND(x)\
    struct timespec x##end;\
    clock_gettime(CLOCK_REALTIME, &x##end); 

#define CTIME(x) double x##time = (x##end.tv_nsec - x##start.tv_nsec);

#define CPRINT(x) printf(#x); printf("(), %f\n", x##time);

#define MSTART(x) auto x##start = std::chrono::steady_clock::now()                                         
#define MEND(x)   auto x##end   = std::chrono::steady_clock::now()                                         
#define MTIME(x)  auto x##time  = std::chrono::duration_cast<std::chrono::nanoseconds>(x##end-x##start).count()
#define MPRINT(x) std::cout << x << "(), " << x##time << std::endl

# Compiler
CXX = g++ 
#CXX=/usr/bin/clang++

# Optimize flags
OPTFLAGS = -O3

# Flags
CPPFLAGS=-I. -I/usr/local/include/eigen3

WARNFLAGS=-Wall -Wno-deprecated

CXXFLAGS = $(OPTFLAGS) $(CPPFLAGS) $(WARNFLAGS) -std=c++11

# Executable source
EXESRCS = main.cc

# Executable object file
EXEOBJS = $(EXESRCS:.cc = .o)

# Executable name
EXEC = darcy

# Sources folder
FOLDER = src/

# External library
GETFEM_PATH = /usr/local/

LIB_PATH =  /usr/lib

LDLIBS = /usr/local/lib/libgetfem.a

LDLIBS += $(GETFEM_LIB) -rdynamic  /usr/lib/libblas.a /usr/lib/liblapack.so.3gf  /usr/lib/x86_64-linux-gnu/libqhull.a

INCLUDE = -I$(GETFEM_PATH) -I$(GETFEM_PATH)/src -I$(GETFEM_PATH)/src/gmm -I./include -I/usr/include 

DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS

# Sources
SRCS = $(wildcard $(FOLDER)*.cc)

# Objects
OBJS = $(SRCS:.cc=.o)

# Headers
HEADERS = $(SRCS:.cc=.h)

# Name file of dependences
DEPEND = make.dep

.PHONY: all clean

all : $(DEPEND) $(OBJS) $(EXEOBJS)
	$(CXX) -g $(OPTFLAGS) -o $(EXEC) $(EXEOBJS) $(OBJS) $(LDLIBS) $(DEF_TAGS) $(INCLUDE) 

$(DEPEND) : $(SRCS) $(EXESRCS)
	$(CXX) -MM $(SRCS) $(EXESRCS) -MF $(DEPEND) 

-include $(DEPEND)

clean :
	-rm $(EXEC) $(OBJS)

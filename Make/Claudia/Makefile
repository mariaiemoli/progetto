# Compiler
CXX ?= g++

# Optimization flags
OPTFLAGS ?= -O3

# Flags
CXXFLAGS = $(OPTFLAGS) -Wall #-std=c++11

# Preprocessor flags
CPPFLAGS += $(OPTFLAGS) -Wall $(INCLUDE) $(DEFINES)

# Executable source
EXESRCS = main.cc

# Executable object file
EXEOBJS = $(EXESRCS:.cc=.o)

# Executable name
EXEC = darcy

# Sources
FOLDER = src/
SRCS = $(wildcard $(FOLDER)*.cc)
OBJS = $(SRCS:.cc=.o)

# External library
GETFEM_PATH ?= /usr/local/include
LIB_PATH ?= /opt/local/lib/

# Linker flags
LDFLAGS += -L$(GETFEM_PATH) -L$(LIB_PATH) 
LIBS = -lblas -llapack -lqhull -lgetfem -latlas

# Includes
INCLUDE = -I$(GETFEM_PATH) -I$(GETFEM_PATH)/getfem -I$(GETFEM_PATH)/gmm -I/opt/local/include/ -I./include

# Defines
DEFINES = -DHAVE_CONFIG_H -DGMM_USES_BLAS

# Name file of dependences
DEPEND = make.dep

.PHONY: all clean distclean

all: $(EXEC)

$(EXEC):$(DEPEND) $(EXEOBJS) $(OBJS)
	$(CXX) -g $(LDFLAGS) $(LIBS) -o $@ $^ 

# the special macros $@ and $^, which are the left and right sides of the :, respectively, to make the overall compilation rule more general.
	
$(DEPEND) : $(SRCS) $(EXESRCS)
	$(CXX) -MM $(SRCS) $(EXESRCS) -MF $(DEPEND) 

-include $(DEPEND)

clean:
	$(RM) $(EXEOBJS) $(OBJS)

distclean: clean
	$(RM) $(EXEC)

# Compiler
CXX = g++

# Optimize flags
OPTFLAGS = -O3

# Flags
CXXFLAGS = $(OPTFLAGS) -Wall

# Executable source
EXESRC = convertMesh.cpp

# Executable name
EXEC = convertMesh

.PHONY: all clean

all :
	$(CXX) $(OPTFLAGS) -o $(EXEC) $(EXESRC)

clean :
	-rm $(EXEC)

EXECUTABLE = SimPLEX

# List of sources to include
include sources.mk

OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
DEPS = $(patsubst %.o, %.d, $(OBJECTS))
MISSING_DEPS = $(filter-out $(wildcard $(DEPS)), $(DEPS))

# Implicit rule for cpp files
#  ‘$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c’ 
# http://www.gnu.org/software/make/manual/html_node/Using-Implicit.html#Using-Implicit

CXX = g++ 

#RM = del 2>nul
# 2>nul suppresses error if file does not exist
# Uncomment for Windows systems

# C PreProcessor flags
CPPFLAGS = -MMD
# Special variable for implicit rule generation
# -MMD generates (-M) the dependency files (*.d) without 
# the system header files (-MM) and does not stop compilation at the 
# preprocessor state (-MMD).

# C++ compiler flags
CXXFLAGS = -g
# Special variable for implicit rule generation
# -g for debugging symbols
 
LDFLAGS = -g
# Special variable for implicit rule generation
# -g for debugging symbols



# This protects against any files called 'all', 'clean', etc. 
.PHONY : all clean rebuild 

all : $(EXECUTABLE)

clean: 
	$(RM) *.o
	$(RM) *.d
	$(RM) $(EXECUTABLE)

rebuild: clean all

# This protects against having a missing dependency file and remakes the object 
ifneq ($(MISSING_DEPS),)
$(MISSING_DEPS) :
	$(RM) $(patsubst %.d, %.o, $@)
endif


$(EXECUTABLE) : $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJECTS)

# The recipes for the dependencies are determined using make's implicit rules
include $(DEPS)


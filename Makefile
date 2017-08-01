# ?= assigns on if the variable is not already defined.
TARGET_EXEC ?= simPLEX

EXEC_DIR ?= ./bin
BUILD_DIR ?= ./build
SRC_DIR ?= ./src

PWD := $(dir $(lastword $(MAKEFILE_LIST)))


ifeq ($(OS),Windows_NT)
	SRCS := $(shell dir $(SRC_DIR)/*.cpp)
	SRCS := $(patsubst $(PWD)%, ./%, $(SRCS))
else
	SRCS := $(shell find $(SRC_DIR) -name *.cpp) #This is specific to linux.
endif
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

ifeq ($(OS),Windows_NT)
	#INC_DIRS := $(shell dir /A:D $(SRC_DIR)/*)
	#INC_DIRS := $(patsubst $(PWD)%, ./%, $(INC_DIRS))
else
	INC_DIRS := $(shell find $(SRC_DIR) -type d)
endif

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

$(EXEC_DIR)/$(TARGET_EXEC): $(OBJS)
	$(MKDIR_P) $(dir $@)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)
	$(RM) -r $(EXEC_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p

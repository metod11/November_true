TARGET = main
SRCS = src/Tools/NiederreiterBaseTwo.cpp src/Tools/StopCondition.cpp src/Methods/bfgs/bfgs.cpp src/Tools/main.cpp src/Methods/powell/powell.cpp src/Methods/dfp/dfp.cpp src/Tools/math.cpp src/Tools/sobolseqgenerator.cpp src/Tools/global_min.cpp src/Methods/nesterov/nesterov.cpp src/Tools/tinyxml2.cpp src/Methods/hessian_free/hessian_free.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"include"

# Поддиректории, в которых также нужно осуществить сборку:
SUBDIRS := test-methods test-NiederreiterBaseTwo

# Если в них нужно осуществлять сборку, то нужно и выполнять make clean
# Преобразуем цели, добавив суффикс, отвечающий за очистку:
CLEAN_SUFFIX := ___clean___
SUBDIRS_CLEAN := $(foreach dir,$(SUBDIRS),$(dir)$(CLEAN_SUFFIX))
EMPTY :=

.PHONY: all $(SUBDIRS) clean

all: $(TARGET) $(SUBDIRS)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
$(SUBDIRS):
		$(MAKE) -C $@
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean: $(SUBDIRS_CLEAN)
		rm -rf $(TARGET) $(OBJS)
$(SUBDIRS_CLEAN):
		$(MAKE) -C $(@:$(CLEAN_SUFFIX)=$(EMPTY)) clean
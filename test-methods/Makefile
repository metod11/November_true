TARGET = main
SRCS = test_methods.cpp ../src/Tools/StopCondition.cpp ../src/Tools/math.cpp ../src/Methods/powell/powell.cpp ../src/Methods/hessian_free/hessian_free.cpp ../src/Methods/nesterov/nesterov.cpp ../src/Methods/bfgs/bfgs.cpp ../src/Methods/dfp/dfp.cpp ../src/Methods/dfp/dfp2.cpp ../src/Methods/AdaMax/AdaMax.cpp ../src/Methods/Adam/Adam.cpp ../src/Methods/Adam/AdamW.cpp ../src/Methods/AMSGrad/AMSGrad.cpp ../src/Methods/Adagrad/Adagrad.cpp ../src/Methods/Adadelta/Adadelta.cpp ../src/Methods/Nadam/Nadam.cpp ../src/Methods/RmsProp/RmsProp.cpp ../src/Methods/powell/powell2.cpp ../src/Methods/bfgs/bfgs2.cpp ../src/Methods/powell/powell21.cpp ../src/Tools/Data.cpp
OBJS = $(SRCS:.cpp=.o)
CXX ?= gcc
CXXFLAGS ?= -pthread -Ofast -std=c++14 -Wall -Wextra -fmax-errors=2 -I"../include"

.PHONY: all clean

all: $(TARGET)
$(TARGET): $(OBJS)
		$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)
 
.c.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -rf $(TARGET) $(OBJS)

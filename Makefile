CC		    = g++
CXXFLAGS    = -std=c++11 -Wall -o3 -flto -g -I/opt/homebrew/opt/boost/include -I/opt/homebrew/Cellar/libomp/16.0.6/include -I/opt/homebrew/Cellar/openblas/0.3.23/include/
LDLIBS		= -lblas -lboost_program_options
LDFLAGS		= -L/opt/homebrew/opt/boost/lib
TARGET		= main
HDRS 		= ShallowWater.h

# Parameters of the 4 test cases specified in handout.
TEST1_PARAMS 	= --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --choice 1
TEST2_PARAMS 	= --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --choice 1
TEST3_PARAMS 	= --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --choice 1
TEST4_PARAMS 	= --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --choice 1

# Setting the "default" target.
default: all
all: $(TARGET)

$(TARGET): $(TARGET).o ShallowWater.o

# A change in the header file requires re-compilation 
# (done as both .cpp include the header file).
%.o: %.cpp $(HDRS)

# These make targets don't create new files but rather run commands.
.PHONY: clean test1 test2 test3 test4

# Cleans working directory.
clean:
	rm -f $(TARGET) *.o


test1: $(TARGET)
	export OMP_NUM_THREADS=1; ./$(TARGET) $(TEST1_PARAMS)

test2: $(TARGET)
	export OMP_NUM_THREADS=1; ./$(TARGET) $(TEST2_PARAMS)
	
test3: $(TARGET)
	export OMP_NUM_THREADS=1; ./$(TARGET) $(TEST3_PARAMS)
	
test4: $(TARGET)
	export OMP_NUM_THREADS=1; ./$(TARGET) $(TEST4_PARAMS)

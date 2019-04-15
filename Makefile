MAKE=make
BIN=penum

CXX = g++ 

SRCDIR=src
CPP_SRCS += \
src/main.cpp 

O_SRCS += \
src/main.o 

OBJS += \
src/main.o 

CPP_DEPS += \
src/main.d 


LANG=-std=c++0x
CXXFLAGS= -O3 -g -Wall -c -fmessage-length=0 -fopenmp -MMD -MP
DEFINES= -DUSE_VECQUEUE
LDFLAGS=-fopenmp -g
INCLUDE=-I../penumlib/src -I../ntl/include
LIBS=-lparntl -lpenumlib -lgmp -lfplll -lm
LIBPATH=-L../penumlib -L../ntl/Release

src/%.o: src/%.cpp
	@echo 'Building file: $<'
	$(CXX) $(LANG) $(DEFINES) $(INCLUDE) $(CXXFLAGS) -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
 
all:
	$(MAKE) --no-print-directory main-build

# Main-build Target
main-build: $(BIN)
	cd penumlib && make
	cd penum && make


clean:
	cd penumlib && make clean
	cd penum && make clean 

echo:

	echo $(OBJS)

-include $(DEPS)

.PHONY: all clean


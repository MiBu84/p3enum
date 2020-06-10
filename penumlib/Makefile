LIB=libpenumlib.a

CXX = g++ 
AR = ar

CPP_SRCS += \
src/Utils.cpp \
src/pEnumeratorDouble.cpp \
src/pruningfunc.cpp 

O_SRCS += \
src/Utils.o \
src/pEnumeratorDouble.o \
src/pruningfunc.o 

OBJS += \
src/Utils.o \
src/pEnumeratorDouble.o \
src/pruningfunc.o 

CPP_DEPS += \
src/Utils.d \
src/pEnumeratorDouble.d \
src/pruningfunc.d 

LANG=-std=c++17
ARFLAGS= -r
CXXFLAGS= -O3 -mavx2 -Wall -g -c -fmessage-length=0 -fopenmp -MMD -MP
DEFINES= -DUSE_VECQUEUE

INCLUDE=-I./include -I../penum/src -I./src

src/%.o: src/%.cpp
	@echo 'Building file: $<'
	$(CXX) $(LANG) $(DEFINES) $(INCLUDE) $(CXXFLAGS) -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'

all: $(LIB)

$(LIB): $(OBJS)
	@echo 'Building target: $@'
	$(AR) $(ARFLAGS) "$(LIB)" $(OBJS)
clean:
	rm -f src/*.d src/*.o $(LIB)

echo:
	echo $(OBJS)

-include $(DEPS)

.PHONY: all clean


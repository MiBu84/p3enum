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
LIBS=-lntl -lpenumlib -lgmp -lfplll -lm
LIBPATH=-L../penumlib -L../ntl/Release

src/%.o: src/%.cpp
	@echo 'Building file: $<'
	$(CXX) $(LANG) $(DEFINES) $(INCLUDE) $(CXXFLAGS) -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
 
all:
	$(MAKE) --no-print-directory main-build

# Main-build Target
main-build: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(LIBPATH) $(LDFLAGS) -o $(BIN) $(OBJS) -Wl,--start-group $(LIBS) -Wl,--end-group
	@echo 'Finished building target: $@'

ntllow:
	cd ../ntl/src/ && ./configure CXXFLAGS="-I../../penumlib/src -g -O3" LDFLAGS="-g" NATIVE=off NTL_THREAD_BOOST=off NTL_THREADS=off NTL_STD_CXX11=off
	cd ../ntl/src/ && make -j 4
	cd ../ntl/src/ && mkdir -p -m 755 ../Release
	cd ../ntl/src/ && cp -p ntl.a ../Release/libparntl.a

reqsincr:
	cd ../ntl/src/ && make -j 4
	cd ../ntl/src/ && mkdir -p -m 755 ../Release
	cd ../ntl/src/ && cp -p ntl.a ../Release/libparntl.a
	cd ../penumlib && make
	rm penum 
	make penum
reqs:
	cd ../ntl/src/ && ./configure CXXFLAGS="-DUSE_VECQUEUE -I../../penumlib/src -g -O3"
	cd ../ntl/src/ && make -j 4
	cd ../ntl/src/ && mkdir -p -m 755 ../Release
	cd ../ntl/src/ && cp -p ntl.a ../Release/libparntl.a
	cd ../penumlib && make
	rm penum
	make

clean:
#	cd ../penumlib && make clean
#	cd ../ntl/src/ && make clean
	cd ../penum && rm -f $(SRCDIR)/*.d $(SRCDIR)/*.o penum

echo:

	echo $(OBJS)

-include $(DEPS)

.PHONY: all clean


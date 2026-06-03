CC = g++
CXXFLAGS += -std=c++11 -fPIC -O2 -fpermissive -I../Common/BaseMath -I../Common/BasePhysics -I../Common/General -I../Common/IDL -I../Common/Multithread -I../Common/Utils -I../sources
# CXXFLAGS += -DDEBUG_CPP11THREAD

LDFLAGS = -shared -lm -pthread

TARGET=../binaries/WWNLFFFReconstruction.so
SOURCES=$(shell echo ../Common/BaseMath/*.cpp ../Common/BasePhysics/*.cpp ../Common/Multithread/*.cpp ../Common/Utils/*.cpp ../sources/*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

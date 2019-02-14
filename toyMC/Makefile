TARGET=resolution
INCS=-I$(shell root-config --cflags)
LIBS=$(shell root-config --libs)

all: $(TARGET)
$(TARGET): $(TARGET).cpp
	g++ -O3 -std=c++17  $(INCS) $(LIBS) -o $(TARGET) $(TARGET).cpp

clean:
		rm -f $(TARGET) core

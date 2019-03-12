TARGET=resolution
INCS=-I$(shell root-config --incdir)
LIBS=$(shell root-config --libs)
all: $(TARGET)
$(TARGET):
	g++ -O3 -std=c++17  $(INCS) $(LIBS) -o $(TARGET) $(TARGET).cpp


TARGET=resolution
INCS=-I$(shell root-config --incdir)
LIBS=$(shell root-config --libs)
all: $(TARGET)
$(TARGET):
	g++ -O3  $(INCS) $(LIBS) -o $(TARGET) $(TARGET).cpp


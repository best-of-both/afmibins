CPP = g++
INCLUDES = 
LIBS = 

LIBDIR = lib/
OBJDIR = obj/
SRCDIR = src/
INCDIR = include/
TESTDIR = test/
BINDIR = bin/

INCLUDES += -I. -I$(INCDIR) -I/opt/local/include/atlas -I/opt/local/include
LIBDIRS += -L/opt/local/lib
LIBS += -lptcblas -latlas -lptlapack -lpthread

OBJECTS = $(patsubst $(SRCDIR)%.cc, $(OBJDIR)%.o, $(shell ls -F $(SRCDIR)*.cc))
TESTS = $(patsubst $(TESTDIR)%.cc, $(BINDIR)test_%, $(shell ls -F $(TESTDIR)*.cc))
TARGET = libatlas++.a

all: $(LIBDIR)$(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)%.o: $(SRCDIR)%.cc $(OBJDIR)
	$(CPP) -g -std=c++11 $(INCLUDES) $(LIBS) -c -o $@ $<

$(LIBDIR)$(TARGET): $(OBJECTS) $(LIBDIR)
	ar cr $(LIBDIR)$(TARGET) $(OBJECTS)

clean:
	rm $(OBJECTS) $(LIBDIR)$(TARGET)

$(BINDIR)test_%: $(TESTDIR)%.cc $(BINDIR)
	$(CPP) -g -std=c++11 $(INCLUDES) $(LIBDIRS) $(LIBS) -L$(LIBDIR) -l$(patsubst lib%.a,%,$(TARGET)) -o $@ $<

test: $(LIBDIR)$(TARGET) $(TESTS)

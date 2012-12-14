#define this module
include ../../MakeDefns

OBJLIB = evalIR.a

all: $(OBJLIB)

# compiling C++ code
%.o: ../src/%.cpp
	$(CXX) $(CXXFLAGS) -I../include -o $@ -c $<

# compiling C code
%.o: ../src/%.c
	$(CC) $(CFLAGS) -I../include -o $@ -c $<

# how to make a library from object files
$(OBJLIB): $(OBJ)
	rm -f $@; $(AR) -rs $@ $(OBJ)

clean:
	rm -f *.o *.a


INCDIR = $(includedir)
ALLHEADER = $(wildcard ../include/*.h*)
install:
	$(INSTALL_DATA) $(ALLHEADER) $(INCDIR)
	$(AR) -qs $(libdir)/$(INDRILIB) $(OBJ)

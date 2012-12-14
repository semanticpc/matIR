include ../MakeDefns
SHARED=
INCPATH=-I./include
LIBPATH=-L./obj
LIBS= $(indriLibrary) obj/evalIR.a -larmadillo
APP1=divEval

all:
	$(MAKE) -C obj -f ../src/Makefile
	$(CXX) $(CXXFLAGS) $(APP1).cpp -o $(APP1) $(LIBPATH) $(LIBS) $(CPPLDFLAGS)

install:

	$(MAKE) install -C obj -f ../src/Makefile
	$(INSTALL_PROGRAM) $(APP1) $(bindir)
clean:
	$(MAKE) clean -C obj -f ../src/Makefile
	rm -f $(APP1)


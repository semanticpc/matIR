include ../MakeDefns
SHARED=
INCPATH=-I../include
LIBPATH=-L../obj
LIBS= $(indriLibrary) ../obj/matIR.a -larmadillo
APP1=docVector
APP2=queryFeatures
APP3=documentFeatures

all:
	$(CXX) $(CXXFLAGS) $(APP1).cpp -o $(APP1) $(LIBPATH) $(LIBS) $(CPPLDFLAGS)
	$(CXX) $(CXXFLAGS) $(APP2).cpp -o $(APP2) $(LIBPATH) $(LIBS) $(CPPLDFLAGS)
	$(CXX) $(CXXFLAGS) $(APP3).cpp -o $(APP3) $(LIBPATH) $(LIBS) $(CPPLDFLAGS)

install:
	$(INSTALL_PROGRAM) $(APP1) $(bindir)
	$(INSTALL_PROGRAM) $(APP2) $(bindir)
	$(INSTALL_PROGRAM) $(APP3) $(bindir)

clean:
	rm -f $(APP1)
	rm -f $(APP2)
	rm -f $(APP3)
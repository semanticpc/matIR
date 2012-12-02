#define this module
MODULE	= matIR

-include ../MakeDefns
OBJLIB = matIR.a

all: $(OBJLIB)


# how to make a library from object files
$(OBJLIB): $(OBJ) 
	rm -f $@; $(AR) -rs $@ $(OBJ)



# compiling C++ code
%.o: ../src/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# compiling C code
%.o: ../src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o *.a

install:
	$(INSTALL_DATA) $(ALLHEADER) $(pkgincludedir)
	$(AR) -rs $(libdir)/$(INDRILIB) $(OBJ)

stamp:
	
# Each source file has a corresponding ".d" file
# that specifies its dependents. 
ifneq (clean, $(MAKECMDGOALS))
ifneq (install, $(MAKECMDGOALS))
ifneq (stamp, $(MAKECMDGOALS))
ifneq ($(strip $(CPPSRC)),)
-include $(patsubst %, ../depend/%, $(CPPSRC:.cpp=.d))
endif

ifneq ($(strip $(CSRC)),)
-include $(patsubst %, ../depend/%, $(CSRC:.c=.d))
endif
endif
endif
endif

# how to generate .d  file from a .C file
../depend/%.d: ../src/%.cpp
	@echo "Generating $@ for $<"
	@set -e; $(CXX) -MM $(CXXFLAGS) $< | sed 's/\($*\)\.o[ :]*/\1.o \.\.\/depend\/\1.d : /g' > $@; [ -s $@ ] || rm -f $@

# how to generate .d  file from a .C file
../depend/%.d: ../src/%.cc
	@echo "Generating $@ for $<"
	@set -e; $(CXX) -MM $(CXXFLAGS) $< | sed 's/\($*\)\.o[ :]*/\1.o \.\.\/depend\/\1.d : /g' > $@; [ -s $@ ] || rm -f $@

# how to generate .d file from a .c file
../depend/%.d: ../src/%.c
	@echo "Generating $@ for $<"
	@set -e; $(CC) -MM $(CFLAGS) $< | sed 's/\($*\)\.o[ :]*/\1.o \.\.\/depend\/\1.d : /g' > $@; [ -s $@ ] || rm -f $@


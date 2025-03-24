
SHELL = /bin/sh

#----------------------

SRC       = $(wildcard *.cpp)
OBJ       = $(SRC:.cpp=.o)
EXE       = Helical_surface
#----------------------

CC        = g++
CCFLAGS   = -g -O 
LINK      = g++
LINKFLAGS = -g -O

DEPFLAGS  = -MM
SIZE      = size
#---------------------
#generate alloy 
#CCFLAGS   += -DALLOY
#CCFLAGS   += -DNT
#CCFLAGS   += -DFCC
#CCFLAGS    += -DSQUARE_LOOP     #create a squre shape dislocation loop 
#CCFLAGS   += -DRESET_DISL_TYPE  #reset element type of dislocator

#----------------------
# LINK TARGET
$(EXE): $(OBJ)
	$(CC) $(LINKFLAGS) -o $@ $(OBJ)
	$(SIZE) $@

# Compilation rules
%.o: %.cpp
	$(CC) $(CCFLAGS) -c $<

%.d: %.cpp
	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@

#individual depdence
DEPDENCE = $(SRC:.cpp=.d)
sinclude $(DEPDENCE)

.PHONY: clean clean-obj
clean-obj:
	rm *.d *.o
clean:
	rm $(EXE) *.d *.o


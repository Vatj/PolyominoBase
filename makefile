#Compiler and Linker
CXX         := g++

#The Target Binary Program
EV_TARGET  := EvolutionSimulator
BR_TARGET   := StochasticAssembler
PR_TARGET   := BulkProcessor


#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := includes
BUILDDIR    := build
TARGETDIR   := bin
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CXXFLAGS    := -std=c++11 -Wall -O3 -fopenmp 
INC         := -I$(INCDIR)
INCDEP      := -I$(INCDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
COMMON_SOURCES     := $(shell find $(SRCDIR) -type f -name graph_*.$(SRCEXT)) $(SRCDIR)/xorshift.$(SRCEXT)
EV_SOURCES   := $(shell find $(SRCDIR) -type f -name evolution_*.$(SRCEXT))
BR_SOURCES := $(shell find $(SRCDIR) -type f -name brute_*.$(SRCEXT)) 
PR_SOURCES := $(shell find $(SRCDIR) -type f -name processing_*.$(SRCEXT)) 

COMMON_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(COMMON_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
EV_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(EV_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
BR_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(BR_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
PR_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(PR_SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Default Make
all: Ev Br Pr

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Pull in dependency info for *existing* .o files
-include $(COMMON_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(EV_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(BR_OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
#$(TARGET): $(OBJECTS)
#	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(TARGET) $^

Ev: $(EV_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(EV_TARGET) $^

Br: $(BR_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(BR_TARGET) $^

Pr: $(PR_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(PR_TARGET) $^

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	@mkdir -p $(dir $(TARGETDIR))
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all clean Ev Br Pr

#Compiler and Linker
CXX         := g++

#The Target Binary Program
EV_TARGET   := EvolutionSimulator
ST_TARGET   := StochasticAssembler
PR_TARGET   := BulkProcessor
PE_TARGET   := ProteinEvolution
GP_TARGET   := GP_Mapping



#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := includes
BUILDDIR    := build
TARGETDIR   := bin
PROFDIR	    := profiling
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CXXFLAGS    := -std=gnu++14 -Wall -Wextra -pedantic -pipe -O3 -fopenmp -march=ivybridge -flto -flto-partition=none #-fprofile-use -fprofile-correction -fprofile-dir=$(PROFDIR)/
INC         := -I$(INCDIR)
INCDEP      := -I$(INCDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
COMMON_SOURCES     := $(shell find $(SRCDIR) -type f -name graph_*.$(SRCEXT)) 
EV_SOURCES   := $(shell find $(SRCDIR) -type f -name evolution_*.$(SRCEXT))
ST_SOURCES := $(shell find $(SRCDIR) -type f -name stochastic_m*.$(SRCEXT)) 
PR_SOURCES := $(shell find $(SRCDIR) -type f -name processing_*.$(SRCEXT))
PE_SOURCES := $(shell find $(SRCDIR) -type f -name interface_*.$(SRCEXT))
GP_SOURCES := $(shell find $(SRCDIR) -type f -name genotype_*.$(SRCEXT))
Core_P_SOURCES := $(shell find $(SRCDIR) -type f -name core_*.$(SRCEXT))


COMMON_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(COMMON_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
EV_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(EV_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
ST_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(ST_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
PR_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(PR_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
PE_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(PE_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
GP_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(GP_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
Core_P_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(Core_P_SOURCES:.$(SRCEXT)=.$(OBJEXT)))


#Default Make
all: Ev St Pr Pe GP

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Pull in dependency info for *existing* .o files
-include $(COMMON_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(EV_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(ST_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(PR_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(PE_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(TEST_OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
#$(TARGET): $(OBJECTS)
#	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(TARGET) $^

Ev: $(EV_OBJECTS) $(COMMON_OBJECTS) $(ST_OBJECTS)
	@mkdir -p $(TARGETDIR)	
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(EV_TARGET) $^

St: $(ST_OBJECTS) $(COMMON_OBJECTS) $(Core_P_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(ST_TARGET) $^

Pr: $(PR_OBJECTS) $(COMMON_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS)  -o $(TARGETDIR)/$(PR_TARGET) $^

Pe: $(PE_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS)  -o $(TARGETDIR)/$(PE_TARGET) $^

GP: $(GP_OBJECTS) $(COMMON_OBJECTS) $(ST_OBJECTS) $(Core_P_OBJECTS)
	@mkdir -p $(TARGETDIR)	
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(GP_TARGET) $^


#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all clean Ev St Pr Pe GP

#Compiler and Linker
CXX         := g++

#The Target Binary Program
EV_TARGET  := EvolutionSimulator
BR_TARGET   := StochasticAssembler
PR_TARGET   := BulkProcessor
PE_TARGET   := ProteinEvolution

TEST_TARGET := TestSuite


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
COMMON_SOURCES     := $(shell find $(SRCDIR) -type f -name graph_*.$(SRCEXT)) $(SRCDIR)/xorshift.$(SRCEXT)
EV_SOURCES   := $(shell find $(SRCDIR) -type f -name evolution_*.$(SRCEXT))
BR_SOURCES := $(shell find $(SRCDIR) -type f -name brute_*.$(SRCEXT)) 
PR_SOURCES := $(shell find $(SRCDIR) -type f -name processing_*.$(SRCEXT))
PE_SOURCES := $(shell find $(SRCDIR) -type f -name interface_*.$(SRCEXT))
TEST_SOURCES := $(shell find $(SRCDIR) -type f -name test_*.$(SRCEXT)) 

COMMON_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(COMMON_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
EV_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(EV_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
BR_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(BR_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
PR_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(PR_SOURCES:.$(SRCEXT)=.$(OBJEXT)))
PE_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(PE_SOURCES:.$(SRCEXT)=.$(OBJEXT)))

TEST_OBJECTS  := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(TEST_SOURCES:.$(SRCEXT)=.$(OBJEXT))) $(BUILDDIR)/xorshift.$(OBJEXT) 


#Default Make
all: Ev Br Pr

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Pull in dependency info for *existing* .o files
-include $(COMMON_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(EV_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(BR_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(PR_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(PE_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(TEST_OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
#$(TARGET): $(OBJECTS)
#	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(TARGET) $^

Ev: $(EV_OBJECTS) $(COMMON_OBJECTS)
	@mkdir -p $(TARGETDIR)	
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(EV_TARGET) $^

Br: $(BR_OBJECTS) $(COMMON_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS) -o $(TARGETDIR)/$(BR_TARGET) $^

Pr: $(PR_OBJECTS) $(COMMON_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS)  -o $(TARGETDIR)/$(PR_TARGET) $^

Pe: $(PE_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS)  -o $(TARGETDIR)/$(PE_TARGET) $^

Test: $(TEST_OBJECTS)
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS)  -o $(TARGETDIR)/$(TEST_TARGET) $^

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
.PHONY: all clean Ev Br Pr Pe Test

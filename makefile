#Compiler and Linker
CC          := g++

#The Target Binary Program
TARGET      := AlphaTree

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      :=
BUILDDIR    := obj
TARGETDIR   := .
RESDIR      := res
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#-Wall -O3
#Flags, Libraries and Includes
CFLAGS      := -lstdc++fs -std=c++11 -fopenmp -Wall -g -O3 #std=gnu9 -pg -ggdb3 -lprofiler -lpthread
LIB         := -fopenmp
INC         := -I /usr/local/include
INCDEP      :=

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
#SOURCES     := src/main.cpp
#OBJECTS     := obj/main.o
SOURCES     := $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: directories $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
#resources: directories
#	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) $(TARGET)

#Pull in dependency info for *existing* .o files
#-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -pg -o $(TARGETDIR)/$(TARGET) $^ $(LIB)


#Compile
obj/main.o: src/main.cpp
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
#	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
#	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
#	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
#	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
#	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp



#Non-File Targets
.PHONY: all remake clean cleaner resources

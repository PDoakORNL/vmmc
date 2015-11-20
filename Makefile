# VMMC Makefile
#
# Author : Lester. O. Hedges
# Email  : lester.hedges+vmmc@gmail.com
# Date   : April 15th 2015

################################# INFO ########################################

# This Makefile can be used to build the VMMC library along with its demos
# and documentation. For detailed information on using the Makefile run make
# without a target, i.e. simply run make at your command prompt.
#
# Makefile style adapted from http://clarkgrubb.com/make-file-style-guide
# Conventions:
#   - Environment and Makefile variables are in upper case, user
#     defined variables are in lower case.
#   - Variables are declared using the immediate assignment operator :=

############################### MACROS ########################################

# the coloring caused me trouble in emacs and isn't good for echo to file either
colorecho = @echo "$2"
boldcolorecho = @echo "$2"
inlinecolorecho = @echo "$2"

############################## VARIABLES ######################################

# Set shell to bash.
SHELL := bash

# Default goal will print the help message.
.DEFAULT_GOAL := help

# Project name.
project := vmmc

# C++ compiler.
CXX := g++-5

# Installation path.
PREFIX := $(HOME)/local

# Path for source files.
src_dir := src

# Path for demo code.
demo_dir := demos

exa_dir := exa

# Path for object files.
obj_dir := obj

# Path for demo object files.
demo_obj_dir := $(demo_dir)/obj

exa_obj_dir := $(exa_dir)/obj

# Path for the library.
lib_dir := lib

# Path for the demonstration library.
demo_lib_dir := $(demo_dir)/lib

exa_lib_dir := $(exa_dir)/lib

# Library header file.
library_header := $(src_dir)/VMMC.h

# Demo header file.
demo_library_header := $(demo_dir)/src/Demo.h

exa_library_header := $(exa_dir)/src/Exa.h

# Generate library target name.
library := $(lib_dir)/lib$(project).a

# Generate demo library target name.
demo_library := $(demo_lib_dir)/libdemo.a

exa_library := $(exa_lib_dir)/libexa.a

# Install command.
install_cmd := install

# Install flags for executables.
iflags_exec := -m 0755

# Install flags for non-executable files.
iflags := -m 0644

outside_includes := -I/home/epd/codes/exafmm/include -I/home/epd/codes/exafmm/vectorclass -I/home/epd/.linuxbrew/include/gsl -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_EXPANSION=4

# Git commit information.
commit := $(shell git describe --abbrev=4 --dirty --always --tags 2> /dev/null)

# Python header file.
# assumes use of locate and one python, unlikely for developers
# python_header := $(shell locate Python.h | grep 2.7 | head -n 1 | awk -F "/Python.h" '{print $1}')
python_header := ${VIRTUAL_ENV}/include/python2.7

# Python library.
#python_library := $(shell locate libpython2.7 | head -n 1)
python_library := /usr/local/Cellar/python/2.7.10/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib

# C++ compiler flags for development build.
cxxflags_devel := -std=c++11 -gdwarf-2 -g3 -Wall -Isrc $(outside_includes) -DCOMMIT=\"$(commit)\" $(OPTFLAGS) -D_GLIBCXX_DEBUG

# C++ compiler flags for release build.
cxxflags_release := -O3 -std=c++11 -funroll-loops -DNDEBUG -Isrc $(outside_includes) -DCOMMIT=\"$(commit)\" $(OPTFLAGS) -fsanitize=address

# Default to release build.
CXXFLAGS := $(cxxflags_release)

# The C++ header, source, object, and dependency files.
vmmc_headers := $(wildcard $(src_dir)/*.h)
vmmc_sources := $(wildcard $(src_dir)/*.cpp)
temp := $(patsubst %.cpp,%.o,$(vmmc_sources))
objects := $(subst $(src_dir),$(obj_dir),$(temp))
-include $(subst .o,.d,$(objects))

# Source files and executable names for demos.
demo_files := $(wildcard $(demo_dir)/*.cpp)
demos := $(patsubst %.cpp,%,$(demo_files))
demo_headers := $(wildcard $(demo_dir)/src/*.h)
demo_headers := $(filter-out $(demo_library_header), $(demo_headers))
demo_sources := $(wildcard $(demo_dir)/src/*.cpp)
temp := $(patsubst %.cpp,%.o,$(demo_sources))
demo_objects := $(subst $(demo_dir)/src,$(demo_obj_dir),$(temp))

# Source files and executable names for exas.
exa_files := $(wildcard $(exa_dir)/*.cpp)
exas := $(patsubst %.cpp,%,$(exa_files))
exa_headers := $(sort $(wildcard $(exa_dir)/src/*.h) $(wildcard $(exa_dir)/src/*.hpp))
exa_headers := $(filter-out $(exa_library_header), $(exa_headers))
exa_sources := $(wildcard $(exa_dir)/src/*.cpp)
temp := $(patsubst %.cpp,%.o,$(exa_sources))
exa_objects := $(subst $(exa_dir)/src,$(exa_obj_dir),$(temp))
ex_exa_libs := /usr/lib/gcc/x86_64-linux-gnu/5/libgomp.so
exa_cxxflags := -DEXAFMM_WITH_OPENMP -fopenmp

# Source files and executable names for Python demos.
python_demo_files := $(wildcard $(demo_dir)/python/*.cpp)
python_demos := $(patsubst %.cpp,%,$(python_demo_files))
python_sources := $(wildcard $(demo_dir)/python/demo/*.py)

# Doxygen files.
dox_files := $(wildcard dox/*.dox)

# ############################### TARGETS #######################################

# Print help message.
.PHONY: help
help:
	$(call boldcolorecho,4,"About")
	@echo " This Makefile can be used to build the $(project) library along with its"
	@echo " associated demos and documentation."
	@echo
	$(call boldcolorecho,4,"Targets")
	@echo " help       -->  print this help message"
	@echo " build      -->  build library and demos (default=release)"
	@echo " devel      -->  build using development compiler flags (debug)"
	@echo " exa    -->  build exafmm_vmmc using debug compiler flags (optimized)"
	@echo " release    -->  build using release compiler flags (optimized)"
	@echo " doc        -->  generate source code documentation with doxygen"
	@echo " clean      -->  remove object and dependency files"
	@echo " clobber    -->  remove all files generated by make"
	@echo " install    -->  install library, demos, and documentation"
	@echo " uninstall  -->  uninstall library, demos, and documentation"
	@echo
	$(call boldcolorecho,4,"Tips")
	@echo " To set a different installation path run"
	@echo "     PREFIX=path make install"
	@echo
	@echo " Additional CXXFLAGS can be passed using OPTFLAGS, e.g."
	@echo "     OPTFLAGS=-Wall make devel"
	@echo
	@echo " To compile an optimised version for pure isotropic systems"
	@echo "     OPTFLAGS=-DISOTROPIC make release"
	@echo
	@echo " Targets can be chained together, e.g."
	@echo "     make release doc"

# Set development compilation flags and build.
devel: CXXFLAGS := $(cxxflags_devel)
devel: build

# Set release compilation flags and build.
release: CXXFLAGS := $(cxxflags_release)
release: build

exa: CXXFLAGS := $(cxxflags_devel) $(exa_cxxflags)
exa: exab

exa-release: CXXFLAGS := $(cxxflags_release) $(exa_cxxflags)
exa-release: exab

# Print compiler flags.
devel release:
	$(call colorecho, 5, "--> CXXFLAGS: $(CXXFLAGS)")

# Save compiler flags to file if they differ from previous build.
# This target ensures that all object files are recompiled if the flags change.
.PHONY: force
.compiler_flags: force
	@echo "$(CXXFLAGS)" | cmp -s - $@ || echo "$(CXXFLAGS)" > $@

# Check that Python header file and library are present.
# This target ensures that all Python demos are recompiled if the .check_python file changes.
.PHONY: force
.check_python: force
	@echo "Python found." | cmp -s - $@ || \
	if [ "$(python_header)" = "" ] || [ "$(python_library)" = "" ] ; then \
        echo "Python not found."; \
	exit 1; \
	else echo "Python found."; \
	fi > $@


# Compile VMMC object files.
# Autodepenencies are handled using a recipe taken from
# http://scottmcpeak.com/autodepend/autodepend.html
$(obj_dir)/%.o: $(src_dir)/%.cpp .compiler_flags
	@echo "--> Building CXX object $*.o"
	$(CXX) $(CXXFLAGS) -c -o $(obj_dir)/$*.o $(src_dir)/$*.cpp
	$(CXX) -MM $(CXXFLAGS) $(src_dir)/$*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$(obj_dir)/$*.o:|' < $*.d.tmp > $(obj_dir)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> $(obj_dir)/$*.d
	@rm -f $*.d.tmp

# Compile demo object files.
# Autodepenencies are handled using a recipe taken from
# http://scottmcpeak.com/autodepend/autodepend.html
$(demo_obj_dir)/%.o: $(demo_dir)/src/%.cpp .compiler_flags
	$(call echo "--> Building CXX object $*.o")
	$(CXX) $(CXXFLAGS) -c -o $(demo_obj_dir)/$*.o $(demo_dir)/src/$*.cpp
	$(CXX) -MM $(CXXFLAGS) $(demo_dir)/src/$*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$(obj_dir)/$*.o:|' < $*.d.tmp > $(demo_obj_dir)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> $(demo_obj_dir)/$*.d
	@rm -f $*.d.tmp

$(exa_obj_dir)/%.o: $(exa_dir)/src/%.cpp .compiler_flags
	@echo "--> Building CXX object $*.o"
	$(CXX) $(CXXFLAGS) -c -o $(exa_obj_dir)/$*.o $(exa_dir)/src/$*.cpp
	$(CXX) -MM $(CXXFLAGS) $(exa_dir)/src/$*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$(obj_dir)/$*.o:|' < $*.d.tmp > $(exa_obj_dir)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> $(exa_obj_dir)/$*.d
	@rm -f $*.d.tmp





# Build the library and demos.
.PHONY: build
build: $(obj_dir) $(demo_obj_dir) $(library) $(demo_library) $(demo_library_header) $(demos) $(python_demos)

# Build the library and demos.
.PHONY: exab
exab: $(obj_dir) $(exa_obj_dir) $(library) $(exa_library) $(exa_library_header) $(exas)


# Create output directory for object and dependency files.
$(obj_dir):
	mkdir $(obj_dir)

# Create output directory for demo object and dependency files.
$(demo_obj_dir):
	mkdir $(demo_obj_dir)

$(exa_obj_dir):
	mkdir $(exa_obj_dir)

# Create demo library header file.
$(demo_library_header): $(demo_headers)
	$(call colorecho,4, "--> Generating CXX library header $(demo_library_header)")
	@echo -e "#ifndef _DEMO_H\n#define _DEMO_H\n" > $(demo_library_header)
	@for i in $(demo_headers);                  \
		do h=`echo $$i | cut -d '/' -f 3`;      \
		echo "#include \"$$h\"";                \
	done | sort -g >> $(demo_library_header)
	@echo -e "\n#endif" >> $(demo_library_header)

# Create exa library header file.
$(exa_library_header): $(exa_headers)
	@echo 4, "--> Generating CXX library header $(exa_library_header)"
	@echo -e "#ifndef _EXA_H\n#define _EXA_H\n" > $(exa_library_header)
	@for i in $(exa_headers);                  \
		do h=`echo $$i | cut -d '/' -f 3`;      \
		echo "#include \"$$h\"";                \
	done | sort -g >> $(exa_library_header)
	@echo -e "\n#endif" >> $(exa_library_header)

# Build the static library.
$(library): $(objects)
	$(call colorecho,1,"--> Linking CXX static library $(library)")
	mkdir -p $(lib_dir)
	ar -rcs $@ $(objects)
	ranlib $@

# Build the static demo library.
$(demo_library): $(demo_objects)
	$(call colorecho, 1, "--> Linking CXX static library $(demo_library)")
	mkdir -p $(demo_lib_dir)
	ar -rcs $@ $(demo_objects)
	ranlib $@

# Build the static demo library.
$(exa_library): $(exa_objects)
	$(call colorecho, 1, "--> Linking CXX static library $(exa_library)")
	mkdir -p $(exa_lib_dir)
	ar -rcsv $@ $(exa_objects)
	ranlib $@


# Compile demonstration code.
$(demos): %: %.cpp $(demo_library_header) $(library) $(demo_library) $(demo_objects)
	$(call colorecho, 1, "--> Linking CXX executable $@")
	-$(CXX) $(CXXFLAGS) -Wfatal-errors -I$(demo_dir)/src $@.cpp $(library) $(demo_library) $(LIBS) $(LDFLAGS) -o $@

$(exas): %: %.cpp  $(exa_library_header) $(library)  $(exa_library) $(exa_objects)
	$(call colorecho, 1, "--> Linking CXX executable $@")
	-$(CXX) $(CXXFLAGS) -Wfatal-errors -I$(exa_dir)/src $@.cpp $(library) $(exa_library) $(LIBS) $(ex_exa_libs) $(LDFLAGS) -o $@


# Compile C++ Python API demonstration code.
$(python_demos): $(python_sources) .check_python .compiler_flags
	$(call colorecho, 1, "--> Linking CXX executable $@")
	-$(CXX) $(CXXFLAGS) -Wfatal-errors -I$(python_header) $@.cpp $(library) $(python_library) $(LIBS) $(LDFLAGS) -o $@

# Build documentation using Doxygen.
doc: $(headers) $(vmmc_headers) $(dox_files)
	$(call colorecho, 4, "--> Generating CXX source documentation with Doxygen")
	doxygen dox/Doxyfile

# Install the library and demos.
.PHONY: install
install: build doc
	$(call colorecho, 3, "--> Installing CXX static library $(library) to $(PREFIX)/lib")
	$(call colorecho, 3, "--> Installing CXX demos $(demos) to $(PREFIX)/share/$(project)-demos")
	$(call colorecho, 3, "--> Installing CXX Doxygen documentation to $(PREFIX)/share/doc/$(project)")
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/lib
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/include/$(project)
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/share/$(project)-demos
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/share/$(project)-demos/python
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/share/$(project)-demos/python/demo
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/share/doc/$(project)
	$(install_cmd) $(iflags) $(library) $(PREFIX)/lib
	$(install_cmd) $(iflags) $(library_header) $(PREFIX)/include/$(project)
	$(install_cmd) $(iflags_exec) $(demos) $(PREFIX)/share/$(project)-demos
	$(install_cmd) $(iflags_exec) $(python_demos) $(PREFIX)/share/$(project)-demos/python
	$(install_cmd) $(iflags) $(python_sources) $(PREFIX)/share/$(project)-demos/python/demo
	cp -r doc/html $(PREFIX)/share/doc/$(project)

# Uninstall the library and demos.
.PHONY: uninstall
uninstall:
	$(call colorecho, 3, "--> Uninstalling CXX static library $(library) from $(PREFIX)/lib")
	$(call colorecho, 3, "--> Uninstalling CXX demos $(demos) from $(PREFIX)/share/$(project)-demos")
	$(call colorecho, 3, "--> Uninstalling CXX Doxygen documentation from $(PREFIX)/share/doc/$(project)")
	rm -f $(PREFIX)/$(library)
	rm -rf $(PREFIX)/include/$(project)
	rm -rf $(PREFIX)/share/$(project)-demos
	rm -rf $(PREFIX)/share/doc/$(project)

# Clean up object and dependecy files.
.PHONY: clean
clean:
	$(call colorecho,6,"--> Cleaning CXX object and dependency files")
	rm -rf $(obj_dir)
	rm -rf $(demo_obj_dir)
	rm -rf $(exa_obj_dir)
	rm -rf $(demo_library_header)
	rm -rf $(lib_dir)
	rm -rf $(demo_lib_dir)
	rm -rf $(exa_lib_dir)
	rm -rf $(exa_library_header)


# Clean up everything produced by make.
.PHONY: clobber
clobber:
	$(call colorecho, 6, "--> Cleaning all output files")
	rm -rf $(obj_dir)
	rm -rf $(demo_obj_dir)
	rm -rf $(lib_dir)
	rm -rf $(demo_lib_dir)
	rm -rf doc
	rm -f $(demos)
	rm -f $(demo_library_header)
	rm -f $(python_demos)
	rm -f $(demo_dir)/python/demo/*.pyc
	rm -f .compiler_flags
	rm -f .check_python

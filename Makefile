# Makefile for p-whap

BIN_DIR=bin/
OBJ_DIR=obj/

FULLWHAP_DEPS=$(OBJ_DIR)fullwhap.o \
	       	$(OBJ_DIR)Matrix.o



INPUTGEN_DEPS=$(OBJ_DIR)inputgen.o

CXX=mpiCC

CXXFLAGS=-O3 -g

INCLUDES=-I.

all: fullwhap inputgen 

${OBJ_DIR}%.o: %.cpp
	@echo '* Compiling $<'
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

.PHONY: fullwhap
fullwhap: $(BIN_DIR)fullwhap

$(BIN_DIR)fullwhap: $(FULLWHAP_DEPS)
	@echo '* Linking $@'
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: inputgen
inputgen: $(BIN_DIR)inputgen

$(BIN_DIR)inputgen: $(INPUTGEN_DEPS)
	@echo '* Linking $@'
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	@echo "Cleaning..."
	rm -rf ${OBJ_DIR} ${BIN_DIR}

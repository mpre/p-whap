# Makefile for p-whap

BIN_DIR=bin/
OBJ_DIR=obj/

WHAP_DEPS=$(OBJ_DIR)whap.o \
	$(OBJ_DIR)Matrix.o \
	$(OBJ_DIR)Bipartition.o

INPUTGEN_DEPS=$(OBJ_DIR)inputgen.o

CXXFLAGS=-O3 -g -fopenmp

INCLUDES=-I.

all: omp-whap inputgen

${OBJ_DIR}%.o: %.cpp
	@echo '* Compiling $<'
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

.PHONY: omp-whap
omp-whap: $(BIN_DIR)omp-whap

$(BIN_DIR)omp-whap: $(WHAP_DEPS)
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
	rm -rf ${OBJ_DIR}* ${BIN_DIR}omp-whap

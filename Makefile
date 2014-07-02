# Makefile for p-whap

BIN_DIR=bin/
OBJ_DIR=obj/

WHAP_DEPS=$(OBJ_DIR)whap.o

all: whap

${OBJ_DIR}%.o: %.cpp
	@echo '* Compiling $<'
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: whap
whap: $(BIN_DIR)whap

$(BIN_DIR)whap: $(WHAP_DEPS)
	@echo '* Linking $@'
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	@echo "Cleaning..."
	rm -rf ${OBJ_DIR} ${BIN_DIR}

# Variables
LIB_DIR = EarthScope/libmseed-3.1.3
LIBMSEED_A = $(LIB_DIR)/libmseed.a
PYX_FILE = libmseed.pyx
C_FILE = libmseed.c
SO_FILE = $(PYX_FILE:.pyx=.so)
BUILD_DIR = build

# Phony targets
.PHONY: default libmseed clean

# Default target
default: libmseed

# Build the Python extension and link with libmseed
libmseed: setup.py $(PYX_FILE) $(LIBMSEED_A)
	@echo "Building Python extension..."
	python3 setup.py build_ext --inplace
	@echo "Cleaning up intermediate files..."
	@rm -f $(C_FILE)
	@rm -Rf $(BUILD_DIR)

# Ensure that libmseed.a is built
$(LIBMSEED_A):
	@echo "Building libmseed..."
	@make -C $(LIB_DIR) libmseed.a

# Clean up the generated files
clean:
	@echo "Cleaning up generated files..."
	@rm -f $(SO_FILE) $(C_FILE)
	@rm -Rf $(BUILD_DIR)
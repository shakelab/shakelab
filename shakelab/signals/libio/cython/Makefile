# Variables
LIB_DIR = EarthScope/libmseed-3.1.3
LIBMSEED_A = $(LIB_DIR)/libmseed.a
PYX_FILE = libmseed.pyx
C_FILE = libmseed.c
SO_FILE = $(PYX_FILE:.pyx=.so)
BUILD_DIR = build

# Command and flag variables
PYTHON = python3
CYTHON = cython
RM = rm -f
RM_RF = rm -Rf
MAKE_CMD = make
CFLAGS = -I$(LIB_DIR)
LDFLAGS = -L$(LIB_DIR)
LDLIBS = -lmseed

# Verbose mode toggle
VERBOSE = 0
Q = $(if $(filter 1,$(VERBOSE)),,@)

# Phony targets
.PHONY: default libmseed clean clean-lib clean-ext install uninstall

# Default target
default: libmseed

# Rule to generate C file from Cython if .pyx file changes
$(C_FILE): $(PYX_FILE)
	$(Q)@echo "Generating C file from Cython..."
	$(Q)$(CYTHON) $(PYX_FILE)

# Build the Python extension and link with libmseed
libmseed: setup.py $(C_FILE) $(LIBMSEED_A)
	$(Q)@echo "Building Python extension..."
	$(Q)$(PYTHON) setup.py build_ext --inplace
	$(Q)@echo "Cleaning up intermediate files..."
	$(Q)$(RM) $(C_FILE)
	$(Q)$(RM_RF) $(BUILD_DIR)

# Ensure that libmseed.a is built with -fPIC
$(LIBMSEED_A):
	$(Q)@echo "Building libmseed with -fPIC..."
	$(Q)@$(MAKE_CMD) -C $(LIB_DIR) CFLAGS="-fPIC $(CFLAGS)" libmseed.a

# Clean up the generated files and library
clean: clean-lib clean-ext

clean-lib:
	$(Q)@echo "Cleaning up libmseed library..."
	$(Q)@$(MAKE_CMD) -C $(LIB_DIR) clean

clean-ext:
	$(Q)@echo "Cleaning up generated files..."
	$(Q)$(RM) $(SO_FILE) $(C_FILE)
	$(Q)$(RM_RF) $(BUILD_DIR)

# Install the Python extension
install:
	$(Q)@echo "Installing Python extension..."
	$(Q)$(PYTHON) setup.py install

# Uninstall the Python extension
uninstall:
	$(Q)@echo "Uninstalling Python extension..."
	$(Q)$(PYTHON) -m pip uninstall -y libmseed


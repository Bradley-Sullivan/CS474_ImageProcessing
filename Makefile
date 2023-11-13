# Compiler and compiler flags
CC = gcc
CFLAGS = -g -O2 -Wall -I./inc -lm
EXENAME = pgm

# Directories
SRC_DIR = src
INCLUDE_DIR = inc
AS_DIR ?= AS1

# List of all parts
PARTS = part1 part2 part3 part4 part5 part6	

# Common source files
COMMON_SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
# Common header files
COMMON_HEADER_FILES := $(wildcard $(INCLUDE_DIR)/*.h)
# Object files for common source files
COMMON_OBJ_FILES := $(patsubst $(SRC_DIR)/%.c,$(SRC_DIR)/%.o,$(COMMON_SRC_FILES))

# Define the default target
all: $(PARTS)

# Rule to compile object files for common source files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(COMMON_HEADER_FILES)
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to build the executable binary for each part
$(PARTS): %: $(AS_DIR)/%/$(EXENAME)

# Rule to build the executable binary for a specific part
$(AS_DIR)/%/$(EXENAME): $(AS_DIR)/%/driver.c $(COMMON_OBJ_FILES) $(COMMON_HEADER_FILES)
	$(CC) $(CFLAGS) $^ -o $@

# Clean rule
clean:
	rm -f $(SRC_DIR)/*.o */*/$(EXENAME)

.PHONY: all clean


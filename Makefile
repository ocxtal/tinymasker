
CC = gcc
GIT = git
RM = rm
MAKE = make
CAT = cat


# install directory
PREFIX = /usr/local


# source directory
SRCDIR = src
TARGET = tinymasker


# compiler flags
OFLAGS = -O3
WFLAGS = -Wall -Wextra -Wshadow
NWFLAGS = $(shell bash -c "if [[ $(CC) = icc* ]]; then echo '-Wno-unused-function'; else echo '-Wno-unused-function -Wno-unused-label -Wno-constant-conversion -Wno-implicit-fallthrough -Wno-missing-field-initializers'; fi")
LDFLAGS = -lpthread $(UTIL_LDFLAGS)
CFLAGS = -std=c99 -march=native $(UTIL_CFLAGS) $(WFLAGS) $(NWFLAGS)
GFLAGS = -g # -fsanitize=memory # -DDEBUG # -fsanitize=address -fsanitize=leak


# intermediate
SRCS = tinymasker.c toml.c re.c

SRCS_INTL = $(addprefix $(SRCDIR)/, $(SRCS))
OBJS_INTL = $(SRCS_INTL:c=o)
DEPS_INTL = $(SRCS_INTL:c=dep)

# default version string is parsed from git tags, otherwise extracted from the source
VERSION = $(shell $(GIT) describe --tags || grep "define TM_VERSION" $(SRCDIR)/tinymasker.c | grep -o '".*"' | sed 's/"//g')


# suffix rule
.c.o:
	$(CC) $(OFLAGS) $(CFLAGS) -c -o $(<:c=o) $<


# rules
all: $(TARGET)

$(TARGET): $(OBJS_INTL)
	$(CC) -o $(TARGET) $(OFLAGS) $(CFLAGS) $(OBJS_INTL) $(LDFLAGS)

debug: $(SRCS_INTL)
	$(CC) -o $(TARGET).debug $(GFLAGS) $(CFLAGS) $(SRCS_INTL) $(LDFLAGS)

clean:
	$(RM) -f $(TARGET) $(TARGET).debug $(OBJS_INTL)

dep: $(DEPS_INTL)
	$(CAT) $^ > $(SRCDIR)/Makefile.dep

$(DEPS_INTL): $(SRCS_INTL)
	$(CC) -MM -MT $(@:dep=o) $(CFLAGS) $(@:dep=c) > $@


# dependencies
-include $(SRCDIR)/Makefile.dep
-include $(SRCDIR)/Makefile.util


CC = gcc
MAKE = make
GIT = git
RM = rm
CP = cp
CAT = cat
MKDIR = mkdir


# install directory
PREFIX = /usr/local


# source directory
SRCDIR = src
TARGET = tinymasker


# compiler flags
OFLAGS   = -O3
WFLAGS   = -Wall -Wextra -Wshadow
NWFLAGS  = $(shell bash -c "if [[ $(CC) = icc* ]]; then echo '-Wno-unused-function -diag-disable=11074,11076'; else echo '-Wno-unused-function -Wno-unused-label -Wno-constant-conversion -Wno-implicit-fallthrough -Wno-missing-field-initializers'; fi")
LDFLAGS  = -lm -lpthread $(UTIL_LDFLAGS)
VFLAGS   = -DTM_VERSION=\"$(VERSION)\" -DTM_COMMIT=\"$(COMMIT)\"
CFLAGS   = -std=c99 -march=native $(UTIL_CFLAGS) $(WFLAGS) $(NWFLAGS) $(VFLAGS)
LTOFLAGS = $(shell bash -c "(printf '\#include \"src/utils/arch.h\"\n\#if defined(_ARCH_GCC_VERSION) && (_ARCH_GCC_VERSION < 470)\n\#error \"\"\n\#endif' | $(CC) $(CFLAGS) -xc - -c -o /dev/null > /dev/null 2>&1) && echo '-flto'")
GFLAGS   = -g
SFLAGS   = $(GFLAGS) -fsanitize=address # -fsanitize=memory -fsanitize=leak


# intermediate
SRCS = dbg.c index.c align.c masker.c patch.c toml.c re.c tinymasker.c

SRCS_INTL = $(addprefix $(SRCDIR)/, $(SRCS))
OBJS_INTL = $(SRCS_INTL:c=o)
DEPS_INTL = $(SRCS_INTL:c=dep)

# default version string is parsed from git tags, otherwise extracted from the source
VERSION = $(shell $(GIT) describe --tags 2>/dev/null || grep "define TM_VERSION" $(SRCDIR)/tinymasker.h | grep -o '".*"' | sed 's/"//g')
COMMIT  = $(shell $(GIT) describe --always --dirty)


# suffix rule
.c.o:
	$(CC) $(OFLAGS) $(LTOFLAGS) $(CFLAGS) -c -o $(<:c=o) $<


# rules
all: $(TARGET)

$(TARGET): $(OBJS_INTL)
	$(CC) -o $(TARGET) $(OFLAGS) $(LTOFLAGS) $(CFLAGS) $(OBJS_INTL) $(LDFLAGS)

debug: $(SRCS_INTL)
	$(CC) -o $(TARGET).d $(GFLAGS) -DDEBUG $(CFLAGS) $(SRCS_INTL) $(LDFLAGS)
	$(CC) -o $(TARGET).g $(GFLAGS) -DDEBUG -DNDEBUG_PRINT -DNDEBUG_BLOCK $(CFLAGS) $(SRCS_INTL) $(LDFLAGS)
	$(CC) -o $(TARGET).s $(SFLAGS) -DDEBUG -DNDEBUG_PRINT -DNDEBUG_BLOCK $(CFLAGS) $(SRCS_INTL) $(LDFLAGS)
	$(CC) -o $(TARGET).t $(OFLAGS) $(CFLAGS) $(SRCS_INTL) -ltcmalloc $(LDFLAGS)

clean:
	$(RM) -f $(TARGET) $(TARGET).g $(TARGET).d $(TARGET).s $(TARGET).t $(OBJS_INTL)

install:
	$(MKDIR) -p $(PREFIX)/bin && $(CP) $(TARGET) $(PREFIX)/bin

dep: $(DEPS_INTL)
	$(CAT) $^ > $(SRCDIR)/Makefile.dep

$(DEPS_INTL): $(SRCS_INTL)
	$(CC) -MM -MT $(@:dep=o) $(CFLAGS) $(@:dep=c) > $@


# dependencies
-include $(SRCDIR)/Makefile.dep
-include $(SRCDIR)/Makefile.util

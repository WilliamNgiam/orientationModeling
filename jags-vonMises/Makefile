# Simple build for the JAGS von Mises module.
#
# Usage:
#   make
#   sudo make install
#
# Override JAGS paths if needed:
#   make JAGS_INCLUDEDIR=/usr/local/include JAGS_LIBDIR=/usr/local/lib

JAGS_INCLUDEDIR ?= /usr/include/JAGS
JAGS_LIBDIR ?= /usr/lib/x86_64-linux-gnu
JAGS_MAJOR ?= 4

CXX ?= g++
CXXFLAGS = -std=c++11 -fPIC -O2 -Wall -I$(JAGS_INCLUDEDIR) -Isrc
LDFLAGS = -shared
LDLIBS = -L$(JAGS_LIBDIR) -ljags -ljrmath

SOURCES = src/vonmises.cc src/distributions/DVonMises.cc
TARGET = vonmises.so
MODDIR = $(JAGS_LIBDIR)/JAGS/modules-$(JAGS_MAJOR)

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

install: $(TARGET)
	install -d $(DESTDIR)$(MODDIR)
	install -m 644 $(TARGET) $(DESTDIR)$(MODDIR)/

clean:
	rm -f $(TARGET)

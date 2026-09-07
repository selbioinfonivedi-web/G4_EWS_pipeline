# Convenience wrapper around src/Makefile, which is upstream's and is left
# unmodified. `make` there already copies Phi, Profile and ppma_2_bmp up into
# this directory.
#
# The build is controlled through CC, CXX and LDFLAGS, e.g.
#
#   make LDFLAGS=-static                                     # static Linux build
#   make CC="clang -arch arm64" CXX="clang++ -arch arm64"     # pick a macOS arch
#
# Do NOT set CFLAGS on the command line: src/Makefile uses `CFLAGS +=`, so a
# command-line value silently discards its -O3 -Wall.

PREFIX ?= /usr/local
BINS = Phi Profile ppma_2_bmp

all:
	$(MAKE) -C src

test: all
	@tests/smoke-test.sh

install: all
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $(BINS) $(DESTDIR)$(PREFIX)/bin/

clean:
	$(MAKE) -C src clean
	rm -f src/ppma_2_bmp/ppma_2_bmp
	rm -f $(BINS)
	rm -f Phi.log Phi.inf.sites Phi.inf.list Profile.log Profile.csv
	rm -f matrix.ppm profile.ppm

.PHONY: all test install clean

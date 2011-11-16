bin_PROGRAMS = fpidump

AM_CFLAGS = -g -Wall -pedantic $(GLIB_CFLAGS)
AM_LDFLAGS = -ltiff $(GLIB_LIBS)

fpidump_SOURCES = fpidump.c

ACLOCAL_AMFLAGS = -I m4

if HAVE_GIT_REPO
ChangeLog: stamp-git
	(cd $(srcdir) && $(GIT) log --pretty=medium --) > $@
stamp-git:
	@last_version=`if test -f stamp-git; then cat stamp-git; else :; fi`; \
	now_version=`cd $(srcdir) && $(GIT) rev-parse HEAD`; \
	if test "X$$now_version" != "X$$last_version"; then \
	  echo $$now_version > stamp-git; \
	  echo "Now at git commit $$now_version..."; \
	fi
.PHONY: stamp-git
endif HAVE_GIT_REPO

DISTCLEANFILES = stamp-git
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.h.in configure \
	build-tools/compile build-tools/depcomp \
	build-tools/install-sh build-tools/missing

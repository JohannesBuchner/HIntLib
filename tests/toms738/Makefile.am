## Process this file with automake to produce Makefile.in

check_PROGRAMS = genin2 genin gfarit gfplys

genin2_SOURCES = genin2.f inlo2.f golo2.f calcc2.f calcv.f setfld.f plymul.f charac.f
genin_SOURCES = genin.f inlo.f golo.f calcc.f calcv.f setfld.f plymul.f charac.f
genin_DEPENDENCIES = gftabs.dat irrtabs.dat
gfarit_SOURCES = gfarit.f charac.f setfld.f itop.f ptoi.f plyadd.f plymul.f plydiv.f
gfplys_SOURCES = gfplys.f charac.f setfld.f itop.f ptoi.f plymul.f
gfplys_DEPENDENCIES = gftabs.dat

CLEANFILES = gftabs.dat irrtabs.dat
EXTRA_DIST = genin2_ori.f genin_ori.f gfplys_ori.f original_source.txt COPYRIGHT.html

gftabs.dat: gfarit@EXEEXT@
	@echo "Creating \"gftabs.dat\""
	test -x ./gfarit@EXEEXT@ || sleep 1
	./gfarit@EXEEXT@

irrtabs.dat: gfplys@EXEEXT@
	@echo "Creating \"irrtabs.dat\""
	test -x ./gfplys@EXEEXT@ || sleep 1
	./gfplys@EXEEXT@



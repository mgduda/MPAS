PWD=$(shell pwd)
EXE_NAME=calc_gwd_fields
NAMELIST_SUFFIX=gwd
override CPPFLAGS += -DCORE_GWD

report_builds:
	@echo "CORE=gwd"

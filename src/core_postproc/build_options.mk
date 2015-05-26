PWD=$(shell pwd)
EXE_NAME=postproc_model
NAMELIST_SUFFIX=postproc
override CPPFLAGS += -DCORE_POSTPROC

report_builds:
	@echo "CORE=postproc"

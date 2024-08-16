TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime

WRAP_PYTHON_TOOL = wrap_python3
WRAP_PYTHON_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PYTHON3_TOOL).sh

SRC_PYTHON = $(wildcard scripts/*.py)

all: bin

bin: $(BIN_PYTHON)
	ln -s -f expression_transform $(TOP_DIR)/bin/expression_transform.py

deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service

#
# Ugh. Calling code wants name with .py.
#
deploy-client: deploy-scripts deploy-libs 
	ln -s -f expression_transform $(TARGET)/bin/expression_transform.py

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs:

include $(TOP_DIR)/tools/Makefile.common.rules

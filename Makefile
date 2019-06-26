TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime


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

deploy-service:

include $(TOP_DIR)/tools/Makefile.common.rules

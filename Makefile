TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime


all: bin

bin: $(BIN_PYTHON)

deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service
deploy-client: deploy-scripts 

deploy-service:

include $(TOP_DIR)/tools/Makefile.common.rules

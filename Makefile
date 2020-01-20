THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))
IMAGE = bayes:1.0
CONTAINER = bayes-container

build: 
	docker build --tag $(IMAGE)  \
	    -f Dockerfile \
        .

run:
	docker run --rm \
		-it \
		-p 7775:22 \
	    -v $(THIS_FILE_DIR):/opt/project \
	    --name $(CONTAINER) \
	    $(IMAGE)

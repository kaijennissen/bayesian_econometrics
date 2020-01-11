#!/bin/bash
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

INPUT_DIR=$BASE_DIR/bayesian_econometric_methods

docker run -it \
 --rm \
 -v $INPUT_DIR:/home/docker \ 
--name bayes \
nuest/mro:latest

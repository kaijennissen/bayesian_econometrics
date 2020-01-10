#!/bin/bash

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


INPUT_DIR=$BASE_DIR/chapter_18/exercise_18_2
echo $BASE_DIRecho $INPUT_DIR

docker run -it \
 --rm \
 -v $INPUT_DIR:/home/docker \ 
--name bayes \
nuest/mro:latest

#!/bin/bash
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

docker run -e PASSWORD=RockerDocker \
  --rm \
  -p 8787:8787 \
  -v $BASE_DIR:/home/rstudio \
  --name bayes \
  bayes:0.1

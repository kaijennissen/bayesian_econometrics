FROM rocker/rstudio:3.6.2

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadbclient-dev \
  libpq-dev \
  libssh2-1-dev \
  unixodbc-dev \
  libsasl2-dev \
  libgsl-dev \
  && install2.r --error \
    --deps TRUE \
    ggplot2 \
    reshape2 \
    Rcpp \
    RcppArmadillo \
    RcppGSL \    
    RcppZiggurat \
    SuppDists

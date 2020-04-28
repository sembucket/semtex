# Stage 1: Build
FROM ubuntu:latest AS build

LABEL description="Container for semtex"

ARG DEBIAN_FRONTEND=noninteractive
RUN TZ="Australia/Melbourne" && ln -snf /usr/share/zoneinfo/${TZ} /etc/localtime && echo ${TZ} > /etc/timezone

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        g++ \
        gcc \
        gfortran \
        bison \
        libopenblas-dev \
        libopenmpi-dev \
        openmpi-bin \
        make \
        cmake

COPY . /usr/src/semtex/

WORKDIR /usr/src/semtex/

RUN mkdir build && \
    cd build && \
    cmake .. && \
    make

CMD [ "/bin/bash" ]

# Stage 2: Deploy smaller image
# FROM alpine:latest AS deploy

# RUN apk --no-cache 

# https://devblogs.microsoft.com/cppblog/using-multi-stage-containers-for-c-development/
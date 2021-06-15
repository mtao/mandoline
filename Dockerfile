# Use an official Python runtime as a parent image
FROM ubuntu:groovy
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y gnupg-utils ca-certificates

# from oneapi containers https://github.com/intel/oneapi-containers/blob/master/images/docker/os-tools-ubuntu18.04/Dockerfile
# add apt repo public key
ARG url=https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
ADD $url /
RUN file=$(basename "$url") && \
    apt-key add "$file" && \
    rm "$file"

# configure the repository
ARG repo=https://apt.repos.intel.com/oneapi
RUN echo "deb $repo all main" > /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update

# Install any needed packages specified in requirements.txt
RUN apt-get install -y git cmake g++ ninja-build
RUN apt-get install -y libgmp3-dev libtbb-dev libboost-thread-dev libmpfr-dev libmpfrc++-dev libeigen3-dev libblas-dev liblapack-dev libcgal-dev protobuf-compiler
RUN apt-get install -y intel-oneapi-tbb-devel
RUN apt-get install -y libsuitesparse-dev

# Set the working directory to /app
WORKDIR /app

RUN mkdir mandoline
#COPY src/ include/ cmake/ examples/ poisson_2d/ wavesim_2d/ fluidsim_2d/ tests/ CMakeLists.txt /app/mandoline/
COPY ./ /app/mandoline/

# Download and compile TriWild
WORKDIR /app/mandoline/build
RUN ls ..
RUN cmake .. -DCMAKE_BUILD_TYPE=Release -DMANDOLINE_USE_OPENGL=OFF -DMTAO_USE_PYTHON=OFF -GNinja
RUN ninja

WORKDIR /data


# docker build -f Dockerfile_cpp . -t skyways
# docker run --rm -it -v $PWD:/app -p 9000:9000 -w /app skyways

FROM ubuntu:24.04

##############################################################################################
# 1) common
##############################################################################################

ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt update && apt install -y apt-utils dialog
ENV DEBIAN_FRONTEND="noninteractive"

##############################################################################################
#
##############################################################################################

RUN apt update && apt install -y build-essential \
                                 libfontconfig1 \
                                 qt5-qmake \
                                 qtbase5-dev \
                                 cmake \
                                 libgtest-dev \
                                 libeigen3-dev \
                                 libflann-dev \
                                 vim \
                                 libcgal-dev \
                                 gdb \
                                 libopenblas-dev \
                                 liblapack-dev \
                                 libarpack2-dev \
                                 libsuperlu-dev \
                                 libboost-all-dev

##############################################################################################
# armadillo
##############################################################################################

COPY armadillo-12.6.6.tar.xz .
RUN tar -xvf armadillo-12.6.6.tar.xz && \
    cd armadillo-12.6.6  && \
    ./configure && \
    make && \
    make install

##############################################################################################
#
##############################################################################################

RUN useradd -m user
USER user

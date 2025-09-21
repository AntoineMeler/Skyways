# docker build -f Dockerfile . -t skyways
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
                                 libboost-all-dev \
                                 libarmadillo-dev \
                                 sudo \
                                 python3 \
                                 python3-pip

##############################################################################################
#
##############################################################################################

RUN pip install numpy --break-system-packages

##############################################################################################
#
##############################################################################################

RUN useradd -m -s /bin/bash user && \
    echo "user ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

USER user

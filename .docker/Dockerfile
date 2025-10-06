# FROM wci-repo.llnl.gov:4567/almalinux:9.5
FROM almalinux:9.5

#RUN apt-get -y install git
#RUN apt-get install python3-pip -y
#RUN apt-get install python3-setuptools
#RUN apt-get install build-essential git git-lfs cmake -y
#RUN pip install -U pip

RUN dnf -y update -x ca-certificates && \
    dnf install -y gcc gcc-c++ cmake git git-lfs openjpeg2 python3-pip wget python3-setuptools 


RUN wget -c http://www.fftw.org/fftw-3.3.10.tar.gz -O - | tar -xz
RUN cd fftw-3.3.10 && ./configure && make && make install && make check
RUN ln -s /usr/lib/aarch64-linux-gnu/libfftw3.so.3 /usr/local/lib/libfftw3.so

COPY ./requirements.txt /requirements.txt
RUN pip install -r requirements.txt
RUN mkdir satist
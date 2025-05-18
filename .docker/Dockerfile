FROM python:3.10-bullseye

RUN apt-get update

RUN apt-get -y install git
RUN apt-get install python3-pip -y
RUN apt-get install python3-setuptools
RUN apt-get install build-essential git git-lfs cmake -y
RUN pip install -U pip

RUN wget -c http://www.fftw.org/fftw-3.3.10.tar.gz -O - | tar -xz
RUN cd fftw-3.3.10 && ./configure && make && make install && make check
RUN ln -s /usr/lib/aarch64-linux-gnu/libfftw3.so.3 /usr/local/lib/libfftw3.so

COPY ./requirements.txt ./requirements.txt
RUN pip install -r requirements.txt
RUN mkdir satist
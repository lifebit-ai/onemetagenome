# Base image
FROM ubuntu:latest

# Maintainer
MAINTAINER Phil Palmer <phil@lifebit.ai>

RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN apt-get update

RUN apt-get install --yes \
 build-essential \
 wget \
 curl

ENV DEBIAN_FRONTEND=noninteractive
RUN wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
RUN tar -xvf KronaTools-2.7.tar
WORKDIR KronaTools-2.7
RUN ./install.pl
#RUN ./updateAccessions.sh
RUN ./updateTaxonomy.sh

WORKDIR /

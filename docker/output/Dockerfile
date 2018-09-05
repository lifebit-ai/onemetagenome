# Base Image
FROM sculpto/python2-and-3:latest

# Metadata
LABEL description="Docker image containing all requirements to convert a csv file to a html table"

# Maintainer
MAINTAINER Phil Palmer <phil@lifebit.ai>

RUN pip install --upgrade csvtotable

RUN mkdir /data

COPY ./output.html /data/
COPY ./Neo_Gene_EC_Map.py /data/
COPY ./swiss_map.tsv /data/

WORKDIR /data/

#CMD ["-h"]

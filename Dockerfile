# Base Image
FROM bioconductor/devel_base2:latest

#Metadata
LABEL description="Docker image containing all requirements for lifebit-ai/onemetagenome phylotree process"

# Maintainer
MAINTAINER Phil Palmer <phil@lifebit.ai>

RUN apt-get -y update && \
	apt-get --yes --force-yes install libcurl4-openssl-dev procps

# R dependencies
RUN Rscript -e "install.packages('taxize')"

RUN mkdir /data && \
	mkdir /data/rscripts/

COPY ./docker_onemetagenome_phylotree.r /data/rscripts/

WORKDIR /data/rscripts/

ENTRYPOINT ["Rscript", "/data/rscripts/docker_onemetagenome_phylotree.r"]

# CMD ["-h"]

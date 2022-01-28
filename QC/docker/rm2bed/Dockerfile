FROM python:3.8-slim

MAINTAINER Julian Lucas, juklucas@ucsc.edu

RUN apt-get update && apt-get install -y git \
	&& pip install --no-cache-dir pandas

## Install rm2bed.py
WORKDIR /opt/

RUN git clone https://github.com/rmhubley/RepeatMasker.git \
	&& cd RepeatMasker \
	&& git reset --hard a58f3130a4fedb7784171a539052277d2cccc690 \
	&& cd .. \
	&& cp RepeatMasker/util/RM2Bed.py . \
	&& chmod +x RM2Bed.py \
	&& rm -rf RepeatMasker

ENV PATH=/opt:$PATH

WORKDIR /data
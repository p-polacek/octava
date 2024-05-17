# syntax=docker/dockerfile:1

FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR .

COPY . .

COPY pip_requirements.txt .
COPY r_requirements.R .

RUN apt-get update && \
    apt-get install -y \
    --no-install-recommends \
    build-essential \
    r-base \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev && \
    pip3 install -r pip_requirements.txt &&  \
    Rscript r_requirements.R 

EXPOSE 5101

ENTRYPOINT FLASK_APP=flask_app.py flask run --host=0.0.0.0

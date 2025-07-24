# syntax=docker/dockerfile:1

FROM openjdk:11-jre-slim

ENV PYTHONUNBUFFERED=1

COPY . .
COPY requirements.txt .
COPY id_rsa /root/.ssh/id_rsa

WORKDIR /octava

# Download external tools
RUN chmod 600 /root/.ssh/id_rsa &&\
    apt-get update &&\
    apt-get install -y git wget tar unzip openssh-client &&\
    ssh-keyscan -H github.com >> /root/.ssh/known_hosts &&\
    git clone git@github.com:p-polacek/octava.git . &&\
    chmod +x start.sh &&\
    apt-get install -y python3.10 python3-pip &&\
    pip install --upgrade pip &&\
    pip install -r ../requirements.txt &&\
    apt-get update &&\
    mkdir data/ &&\
    mkdir tmp/ &&\
    mkdir ./modules/tools &&\
    cd ./modules/tools/ &&\
    wget https://cnsgenomics.com/software/gctb/download/gctb_2.5.2_Linux.zip -O gctb.zip &&\
    unzip gctb.zip &&\
    wget https://github.com/rgcgithub/regenie/releases/download/v4.1/regenie_v4.1.gz_x86_64_Linux.zip -O regenie.zip &&\
    unzip regenie.zip &&\
    mv regenie_v4.1.gz_x86_64_Linux regenie &&\
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip -O plink1.9.zip &&\
    unzip plink1.9.zip &&\
    wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241124.zip -O plink2.zip &&\
    unzip plink2.zip &&\
    wget http://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar -O beagle.jar &&\
    git clone https://github.com/Menggg/BLINK.git &&\
    unzip BLINK/blink_linux.zip -d BLINK/ &&\
    apt-get update &&\
    apt-get install -y emboss gemma ncbi-blast+ cd-hit

#    wget https://www.haplotype.org/download/pimass-latest.tar.gz -O pimass.tar.gz &&\
#    tar -zxf pimass.tar.gz &&\
#    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz -O blast.tar.gz &&\
#    mkdir blast/ &&\
#    tar -zxf blast.tar.gz -C blast/ --strip-components=1 &&\
#    rm blast.tar.gz &&\
#    wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz -O cdhit.tar.gz &&\
#    wget ftp://emboss.open-bio.org:21/pub/EMBOSS/EMBOSS-6.6.0.tar.gz -O emboss.tar.gz &&\
#    mkdir emboss/ &&\
#    tar -zxf emboss.tar.gz -C emboss/ --strip-components=1 &&\
#    rm emboss.tar.gz &&\
#    mkdir cdhit/ &&\
#    tar -zxf cdhit.tar.gz -C cdhit/ --strip-components=1 &&\
#    cd cdhit/ &&\
#    make &&\
#    rm cdhit.tar.gz &&\
#    wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz -O gemma.gz &&\
#    gunzip gemma.gz &&\

EXPOSE 5101

# ENTRYPOINT FLASK_APP=flask_app.py flask run --host=0.0.0.0 --port 5101

ENTRYPOINT ["/bin/bash", "./start.sh"]
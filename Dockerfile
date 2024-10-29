# syntax=docker/dockerfile:1

FROM python:3.10-slim-bookworm

ENV PYTHONUNBUFFERED=1

WORKDIR .

COPY . .

COPY requirements.txt .

# Download external tools
RUN apt-get update &&\
    apt-get install -y openjdk-11-jre emboss ncbi-blast+ cd-hit plink1.9 plink2 beagle gemma &&\
    pip install -r requirements.txt &&\
    cd ./modules/tools/ &&\
    wget https://www.haplotype.org/download/pimass-latest.tar.gz -O pimass.tar.gz &&\
    tar -zxf pimass.tar.gz &&\
    git clone https://github.com/Menggg/BLINK.git &&\
    unzip BLINK/bin.zip -d BLINK/


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
#    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip -O plink1.zip &&\
#    unzip plink1.zip &&\
#    wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20241020.zip -O plink2.zip &&\
#    unzip plink2.zip &&\
#    wget http://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar -O beagle.jar &&\
#    wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz -O gemma.gz &&\
#    gunzip gemma.gz &&\
 
EXPOSE 5101

ENTRYPOINT FLASK_APP=flask_app.py flask run --host=0.0.0.0

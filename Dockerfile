# syntax=docker/dockerfile:1

##########################
# Stage 1: Builder
##########################
FROM openjdk:11-jre-slim AS builder

COPY .id /root/.ssh/.id
COPY requirements.txt .

WORKDIR /octava

RUN chmod 600 /root/.ssh/.id &&\
    apt-get update &&\
    apt-get install -y git wget tar unzip openssh-client &&\
    apt-get install -y python3.10 python3-pip &&\
    pip install --upgrade pip &&\
    pip install -r ../requirements.txt &&\
    apt-get update &&\
    ssh-keyscan -H github.com >> /root/.ssh/known_hosts &&\
    echo "Host github.com\n  IdentityFile /root/.ssh/.id\n  StrictHostKeyChecking yes" > /root/.ssh/config &&\
    git clone git@github.com:p-polacek/octava.git . &&\
    chmod +x start.sh &&\
    cd modules/tools/ &&\
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz -O blast.tar.gz &&\
    mkdir blast/ &&\
    tar -zxf blast.tar.gz -C blast/ --strip-components=1  ncbi-blast-2.16.0+/bin/blastp ncbi-blast-2.16.0+/bin/blastdbcmd &&\
    wget https://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip -O snpeff.zip &&\
    unzip snpeff.zip snpEff/snpEff.jar &&\
    wget https://cnsgenomics.com/software/gctb/download/gctb_2.5.2_Linux.zip -O gctb.zip &&\
    unzip -j gctb.zip gctb_2.5.2_Linux/gctb -d . &&\
    wget https://github.com/rgcgithub/regenie/releases/download/v4.1/regenie_v4.1.gz_x86_64_Linux.zip -O regenie.zip &&\
    unzip regenie.zip &&\
    mv regenie_v4.1.gz_x86_64_Linux regenie &&\
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip -O plink1.9.zip &&\
    unzip plink1.9.zip plink &&\
    wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250806.zip -O plink2.zip &&\
    unzip plink2.zip plink2 &&\
    wget http://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar -O beagle.jar &&\
    rm blast.tar.gz gctb.zip regenie.zip plink1.9.zip plink2.zip snpeff.zip &&\
    git clone https://github.com/Menggg/BLINK.git &&\
    unzip BLINK/blink_linux.zip -d BLINK/ &&\
    mv BLINK/blink_linux ./blink &&\
    rm -rf BLINK/ &&\
    apt-get update &&\
    # apt-get install -y gemma cd-hit &&\
    apt-get purge -y gcc g++ make manpages && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /usr/share/doc/* /usr/share/man/* /root/.cache/pip
    
##########################
# Stage 2: Runtime
##########################

FROM openjdk:11-jre-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /octava

COPY --from=builder /usr/local /usr/local
COPY --from=builder /octava /octava

COPY .id /root/.ssh/.id

RUN chmod 600 /root/.ssh/.id &&\
    apt-get update && apt-get install -y \
    git openssh-client python3.10 python3-pip gemma cd-hit &&\
    ssh-keyscan -H github.com >> /root/.ssh/known_hosts &&\
    echo "Host github.com\n  IdentityFile /root/.ssh/.id\n  StrictHostKeyChecking yes" > /root/.ssh/config &&\
    chmod +x start.sh process_external.py && mkdir -p data tmp

EXPOSE 5101
ENTRYPOINT ["/bin/bash", "./start.sh"]

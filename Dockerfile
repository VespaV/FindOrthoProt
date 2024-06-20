# Используем базовый образ Ubuntu 24.04
FROM ubuntu:24.04

# Обновление и установка базовых пакетов
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        unzip \
        git \
        build-essential \
        cmake \
        flex \
        bison \
        libgmp3-dev \
        python3.8 \
        python3-pip \
        perl \
        zlib1g \
        ghostscript \
        tcsh \
        make \
        gcc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Установка RAxML-NG
RUN wget https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip && \
    unzip raxml-ng_v1.2.2_linux_x86_64.zip && \
    mv raxml-ng /usr/local/bin && \
    rm raxml-ng_v1.2.2_linux_x86_64.zip

# Установка MEME Suite
WORKDIR /usr/src/app
COPY . .

RUN wget https://meme-suite.org/meme/meme-software/5.5.5/meme-5.5.5.tar.gz && \
    tar zxf meme-5.5.5.tar.gz && \
    cd meme-5.5.5 && \
    ./configure --prefix=/root/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && \
    make install

# Установка дополнительных системных пакетов
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        muscle3 \
        hmmer \
        fasttree \
        ncbi-blast+ \
        clustalo \
        mafft \
        libxml-parser-perl \
    && rm -rf /var/lib/apt/lists/*

# Установка BLAST
RUN mkdir -p /usr/src/BLAST && \
    cd /usr/src/BLAST && \
    wget -c "ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" && \
    gunzip uniprot_sprot.fasta.gz && \
    makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot

# Установка зависимостей Python
COPY requirements.txt .
RUN python3 -m pip install -r requirements.txt

# Экспозиция порта и установка переменных окружения
EXPOSE 8000
ENV DEBUG=1 \
    SECRET_KEY="0930d30j9jd09j09j109fj01j9f" \
    ALLOWED_HOSTS="localhost,127.0.0.1:80,95.163.223.58,https://findorthoprot.ru,http://findorthoprot.ru,findorthoprot.ru" \
    PATH="/root/meme/bin:/root/meme/libexec/meme-5.5.5:$PATH"

# Запуск приложения
CMD ["python3", "manage.py", "runserver", "0.0.0.0:80"]

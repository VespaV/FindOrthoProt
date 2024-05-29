FROM ubuntu:latest

# Обновление пакетного менеджера и установка необходимых пакетов
RUN apt-get update && \
    apt-get install -y wget python3.8 python3-pip perl zlib1g ghostscript bash tcsh make gcc && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Установка MEME Suite
WORKDIR /usr/src/app
# Копирование файлов приложения
COPY . .

RUN wget https://meme-suite.org/meme/meme-software/5.5.5/meme-5.5.5.tar.gz && \
    tar zxf meme-5.5.5.tar.gz && \
    cd meme-5.5.5 && \
    ./configure --prefix=/root/meme --enable-build-libxml2 --enable-build-libxslt && \
    make && \
    make install

# Установка дополнительных системных пакетов
RUN apt-get update && \
    apt-get install -y muscle3 hmmer fasttree ncbi-blast+ clustalo mafft raxml

RUN cd /usr/src/ && mkdir BLAST && cd BLAST && wget -c "ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" && \
    gunzip uniprot_sprot.fasta.gz && \
    makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot

# Добавление пути MEME Suite в переменную PATH
#ENV PATH="/usr/local/meme/bin:${PATH}"

# Установка зависимостей Python
COPY requirements.txt .

RUN python3 -m pip install -r requirements.txt --break-system-packages

RUN apt-get install -y libxml-parser-perl
# Экспозиция портов и установка переменных окружения
EXPOSE 8000
ENV DEBUG=1
ENV SECRET_KEY="0930d30j9jd09j09j109fj01j9f"
ENV ALLOWED_HOSTS="localhost, 127.0.0.1:80, 95.163.223.58, findorthoprot.ru"

ENV PATH="/root/meme/bin:/root/meme/libexec/meme-5.5.5:$PATH"

# Установка переменной PATH и запуск приложения
CMD ["python3", "manage.py", "runserver", "0.0.0.0:80"]
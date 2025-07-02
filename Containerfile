FROM ubuntu:21.04
RUN ln -snf /usr/share/zoneinfo/America/Los_Angeles /etc/localtime && echo America/Los_Angeles > /etc/timezone && \
    apt-get update && \
    apt-get install --no-install-recommends -y ca-certificates && \
    echo 'deb [trusted=yes] https://infra-east-1.s3.amazonaws.com/repositories/apt focal main' > /etc/apt/sources.list.d/atomic-ai.list && \
    apt-get update && \
    apt-get install --no-install-recommends -y awscli bowtie2 fastqc gnupg python3-minimal python3-pip ngmerge samtools umicollapse && \
    rm -rf /var/lib/apt/lists/* && \
    gpg --list-keys # Creates some directories that Jill needs
RUN pip3 install --upgrade jill pip && \
    jill install --confirm
COPY requirements.jl /app/
COPY requirements.txt /app/
RUN pip3 install -r /app/requirements.txt && \
    julia /app/requirements.jl
COPY * /app/
CMD python3 /app/run_with_s3.py

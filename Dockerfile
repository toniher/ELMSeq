FROM biocorecrg/debian-perlbrew-pyenv

VOLUME /input
VOLUME /output

# Locales, for printing
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base r-base-dev

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

RUN mkdir -p /usr/local/bin
COPY ELMSeq.py /usr/local/bin/ELMSeq.py

# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*



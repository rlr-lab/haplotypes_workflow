# -----------
# Stage 1
# -----------

FROM ubuntu:22.04 AS builder

## ----- System dependencies -----
RUN apt-get update
RUN apt-get install -y --no-install-recommends\
    build-essential \
    wget \
    git \
    curl \
    python3 python3-pip python3-dev python3-venv python3-virtualenv python3-all-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    pkg-config \
    autoconf automake make cmake gcc libtool \
    perl libcurl4-gnutls-dev libssl-dev libdeflate-dev \
    # samtools
    samtools \
    && apt-get clean

## ----- minimap2 -----
RUN git clone https://github.com/lh3/minimap2.git /opt/minimap2 && cd /opt/minimap2 && make -j4

## ----- seqtk -----
RUN git clone https://github.com/lh3/seqtk.git /opt/seqtk && cd /opt/seqtk && make -j4

## ----- htslib -----
RUN mkdir -p /opt/htslib && wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2 -O /opt/htslib/htslib-1.22.1.tar.bz2 && cd /opt/htslib && tar -xvf htslib-1.22.1.tar.bz2 && cd htslib-1.22.1/ && ./configure && make -j4 && make install

## ----- ivar -----
RUN git clone https://github.com/andersen-lab/ivar.git /opt/ivar && cd /opt/ivar && ./autogen.sh && ./configure && make -j4 && make install

## ----- fastplong -----
RUN wget http://opengene.org/fastplong/fastplong -O /usr/local/bin/fastplong && chmod +x /usr/local/bin/fastplong

## ----- Flye -----
RUN git clone https://github.com/fenderglass/Flye /opt/flye && cd /opt/flye && make

# -----------
# Stage 2
# -----------
FROM ontresearch/medaka

LABEL org.opencontainers.image.title="align-tools"
LABEL org.opencontainers.image.description="Long-read viral alignment toolkit"
LABEL org.opencontainers.image.version="1.0"
LABEL org.opencontainers.image.authors="Seth Borrowman"

## ----- System dependencies -----
RUN apt-get update
RUN apt-get install -y --no-install-recommends\
    python3 python3-pip python3-dev \
    perl \
    # samtools
    samtools \
    zlib1g libbz2-1.0 liblzma5 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## ----- minimap2 -----
COPY --from=builder /opt/minimap2/minimap2 /usr/local/bin/

## ----- seqtk -----
COPY --from=builder /opt/seqtk/seqtk /usr/local/bin/

## ----- ivar -----
COPY --from=builder /usr/local/bin/ivar /usr/local/bin/

## ----- fastplong -----
COPY --from=builder /usr/local/bin/fastplong /usr/local/bin/

## ----- Flye -----
COPY --from=builder /opt/flye/flye /usr/local/bin

## ----- Strip binaries -----
RUN strip /usr/local/bin/* || true

ENTRYPOINT ["/bin/bash"]

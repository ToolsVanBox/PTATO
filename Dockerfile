FROM ubuntu:22.04

WORKDIR /ptato

COPY Dockerfile .
RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone

RUN apt-get update && apt-get install -y \
   software-properties-common  && \
   add-apt-repository -y ppa:apptainer/ppa && \
   apt-get -y update && \
   apt-get install -y apptainer

RUN apt-get update && apt-get install -y \
   software-properties-common \
   build-essential \
   libssl-dev \
   uuid-dev \
   libgpgme11-dev \
   squashfs-tools \
   libseccomp-dev \
   pkg-config \
   wget \
   curl \
   openjdk-18-jdk \
   git tabix make \
   libncurses-dev libbz2-dev liblzma-dev r-base r-base-core r-recommended r-base-dev \
   r-cran-rgl r-cran-rjags r-cran-snow r-cran-ggplot2 r-cran-igraph \
   r-cran-lme4 r-cran-rjava r-cran-devtools r-cran-roxygen2 r-cran-rjava \
   libxml2-dev libcairo2-dev libsqlite3-dev libmariadbd-dev \
   libpq-dev libssh2-1-dev unixodbc-dev libcurl4-openssl-dev libssl-dev coinor-libcbc-dev coinor-libclp-dev libglpk-dev

RUN  Rscript -e 'install.packages("BiocManager")'
RUN  Rscript -e 'install.packages("randomForest")'
RUN  Rscript -e 'BiocManager::install(c("BSgenome", "copynumber", "MutationalPatterns", "StructuralVariantAnnotation", "VariantAnnotation", "BSgenome.Hsapiens.UCSC.hg38"))'

RUN git clone https://github.com/nextflow-io/nextflow
RUN cd nextflow && make compile && make pack && make install
WORKDIR /ptato/nextflow
RUN echo $PWD
RUN chmod 777 ./nextflow

CMD bash
# Note: this is the last base version supporting ubuntu focal, not jammy
FROM rocker/rstudio:4.2.1

ARG GH_PAT='NOT_SET'

##  Add Bioconductor system dependencies
RUN wget -O install_bioc_sysdeps.sh https://raw.githubusercontent.com/Bioconductor/bioconductor_docker/master/bioc_scripts/install_bioc_sysdeps.sh \
    && bash ./install_bioc_sysdeps.sh \
    && rm ./install_bioc_sysdeps.sh

# NOTE: if anything breaks the dockerhub build cache, you will probably need to build locally and push to dockerhub.
# After the cache is in place, builds from github commits should be fast.
# NOTE: locales / locales-all added due to errors with install_deps() and special characters in the DESCRIPTION file for niaid/dsb \
RUN apt-get update -y \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		libhdf5-dev \
		libpython3-dev \
		python3-pip \
        locales \
        locales-all \
	&& pip3 install umap-learn demuxEM scikit-learn\
    && pip3 install git+https://github.com/bbimber/GMM-Demux \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

ENV RETICULATE_PYTHON=/usr/bin/python3
ENV USE_GMMDEMUX_SEED=1

# NOTE: this is required when running as non-root. Setting MPLCONFIGDIR removes a similar warning.
ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp

# Let this run for the purpose of installing/caching dependencies
RUN if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
    && Rscript -e "install.packages(c('devtools', 'stringi', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)" \
	&& echo "local({options(repos = BiocManager::repositories())})" >> ~/.Rprofile \
	&& echo "Sys.setenv(R_BIOC_VERSION=as.character(BiocManager::version()));" >> ~/.Rprofile \
    # NOTE: this was added to avoid the build dying if this downloads a binary built on a later R version
    && echo "Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true');" >> ~/.Rprofile \
    # To avoid pthread_create() error. See: https://github.com/bmbolstad/preprocessCore/issues/1 and https://github.com/bmbolstad/preprocessCore/issues/12
    && Rscript -e "remotes::install_github('bmbolstad/preprocessCore', dependencies = T, upgrade = 'always', configure.args = '--disable-threading')" \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# This should not be cached if the files change
ADD . /cellhashR

#NOTE: do manual install of fixed DelayedMatrixStats until new version is pushed
RUN cd /cellhashR \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
    # NOTE: currently BioConductor reports preprocessCore version 1.56.0, while the github repo reports 1.55.2. Therefore do this manual install after installing cellhashR:
    && Rscript -e "remotes::install_github('bmbolstad/preprocessCore', dependencies = T, upgrade = 'always', force = TRUE, configure.args = '--disable-threading')" \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds

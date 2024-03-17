FROM ghcr.io/bimberlabinternal/discvr-base:latest

ARG GH_PAT='NOT_SET'

# This should not be cached if the files change
ADD . /cellhashR

#NOTE: do manual install of fixed DelayedMatrixStats until new version is pushed
RUN cd /cellhashR \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
    # To avoid pthread_create() error. See: https://github.com/bmbolstad/preprocessCore/issues/1 and https://github.com/bmbolstad/preprocessCore/issues/12
    && Rscript -e "remotes::install_github('bmbolstad/preprocessCore', dependencies = T, upgrade = 'always', force = TRUE, configure.args = '--disable-threading')" \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds

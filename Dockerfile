FROM continuumio/miniconda3:latest

LABEL maintainer="Amina Lemsara <lemsaraamina@gmail.com>"

# Set up locales properly
RUN apt-get update && \
    apt-get install --yes --no-install-recommends locales && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Use bash as default shell, rather than sh
ENV SHELL /bin/bash

WORKDIR /tmp

# download the conda environment OR copy from local directory
COPY environment.yml /tmp/environment.yml

# install conda environment
RUN conda update conda --yes && \
    conda env update -v -n root --file /tmp/environment.yml && \
#     conda list --name root && \
    conda clean --all --yes && \
    conda clean -tipy && \
    rm /tmp/environment.yml

WORKDIR /ds
COPY src ./src
RUN chmod +x ./src

# clear tmp if there is anything in there...
RUN rm -rf /tmp/*

# # set bash as the default entry point
ENTRYPOINT [ "/bin/bash", "-c" ]

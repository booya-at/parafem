FROM condaforge/linux-anvil:latest

# install some libs
RUN yum install -y mesa-libGL-devel

ARG UID=1000
COPY entrypoint.sh /opt/docker/bin/entrypoint
#RUN unalias cp

ENV PATH="/opt/conda/bin:${PATH}"
ENV HOST_USER_ID=${UID}
FROM ubuntu:20.04

USER root
ENV DEBIAN_FRONTEND=noninteractive 
ENV TZ=Europe/Prague

RUN apt update
RUN apt install -y python3-rdkit 
RUN apt install -y python3-notebook
RUN apt install -y python3-pil

ENV HOME=/work

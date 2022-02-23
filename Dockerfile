FROM ubuntu:20.04

USER root
ENV DEBIAN_FRONTEND=noninteractive 
ENV TZ=Europe/Prague

RUN apt update
RUN apt install -y python3-rdkit 
RUN apt install -y python3-notebook
RUN apt install -y python3-pil

RUN apt install -y python3-pip
RUN pip3 install py3dmol

RUN apt install -y git
RUN mkdir -p /usr/local/share/jupyter/nbextensions && cd /usr/local/share/jupyter/nbextensions && git clone https://github.com/lambdalisue/jupyter-vim-binding vim_binding && jupyter nbextension enable vim_binding/vim_binding --sys-prefix

ENV HOME=/work

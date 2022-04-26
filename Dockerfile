FROM quay.io/pypa/manylinux2014_x86_64:latest

RUN curl -sSf https://sh.rustup.rs | sh -s -- -y

ENV PATH "$HOME/.cargo/bin:$PATH"


#   docker build -t manylinux-rs:latest .
#   docker run -i -v `pwd`:/io -t manylinux-rs:latest /bin/bash

# https://github.com/k3yavi/vpolo
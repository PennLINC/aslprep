# ASLPrep Docker Container Image distribution
#
# Main image: FROM aslprep_build + install aslprep package.
# Base image is built from Dockerfile.base (like fMRIPrep).
#
# MIT License
#
# Copyright (c) The NiPreps Developers
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

ARG BASE_IMAGE=pennlinc/aslprep_build:0.0.20
FROM ${BASE_IMAGE}

# Install aslprep
COPY . /src/aslprep

ARG VERSION=0.0.1

RUN echo "${VERSION}" > /src/aslprep/aslprep/VERSION && \
    echo "include aslprep/VERSION" >> /src/aslprep/MANIFEST.in && \
    pip install --no-cache-dir "/src/aslprep[doc,maint,test]"

RUN find $HOME -type d -exec chmod go=u {} + && \
    find $HOME -type f -exec chmod go=u {} + && \
    rm -rf $HOME/.npm $HOME/.conda $HOME/.empty

RUN ldconfig
WORKDIR /tmp/

ENTRYPOINT ["/opt/conda/envs/aslprep/bin/aslprep"]

ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="aslprep" \
      org.label-schema.description="ASLPrep: A Robust Preprocessing Pipeline for ASL Data" \
      org.label-schema.url="https://aslprep.readthedocs.io/" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/PennLINC/aslprep" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"

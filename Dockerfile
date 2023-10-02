FROM pennlinc/aslprep_build:0.0.2

# Install aslprep
COPY . /src/aslprep

ARG VERSION=0.0.1

# Force static versioning within container
RUN echo "${VERSION}" > /src/aslprep/aslprep/VERSION && \
    echo "include aslprep/VERSION" >> /src/aslprep/MANIFEST.in && \
    pip install --no-cache-dir "/src/aslprep[doc,maint,test]"

RUN find $HOME -type d -exec chmod go=u {} + && \
    find $HOME -type f -exec chmod go=u {} + && \
    rm -rf $HOME/.npm $HOME/.conda $HOME/.empty

RUN ldconfig
WORKDIR /tmp/

ENTRYPOINT ["/usr/local/miniconda/bin/aslprep"]

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

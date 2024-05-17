FROM --platform=linux/amd64 python:3.9
RUN pip install scikit-learn==0.24.2 biopython rich rich-argparse
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz
RUN tar -xvf hmmer.tar.gz
RUN cd hmmer-3.4 && ./configure --prefix=/usr/local && make && make install
RUN rm -rf hmmer.tar.gz hmmer-3.4
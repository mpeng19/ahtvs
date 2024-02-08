FROM ubuntu:latest
RUN apt-get update --fix-missing && apt-get install -yq --no-install-recommends \
    wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git && apt-get clean && \
    rm -rf /var/lib/apt/lists/*


    
# download conda
#RUN printf '. /opt/conda/etc/profile.d/conda.sh\nconda activate' > /etc/profile.d/conda.sh && \
#    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
#    rm ~/miniconda.sh

RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN printf '. /opt/conda/etc/profile.d/conda.sh\nconda activate' > /etc/profile.d/conda.sh

RUN /opt/conda/bin/conda update conda --yes


ADD requirements.txt requirements.txt
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate && \
    pip install --upgrade pip && \
    conda install --yes -c rdkit rdkit python=3 && \
    conda install --yes psycopg2 matplotlib pandas jupyter
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate && \
    conda install --yes -c conda-forge --no-chan-pri django django-debug-toolbar \
    django-extensions django-filter django-guardian djangorestframework && \
    pip install -r requirements.txt
#    pip install -e /a2g2/djangochem

# setup environment
#RUN pip install --upgrade pip
#ADD requirements.txt requirements.txt

#RUN conda update conda --yes
#RUN conda install --yes conda-build  # required just for the conda develop call below
#RUN conda install --yes --channel https://conda.anaconda.org/rdkit rdkit==2016.09.4 python=3.5
#RUN conda install --yes psycopg2 matplotlib pandas jupyter
#RUN conda install --yes -c conda-forge django django-debug-toolbar \
#    django-extensions django-filter django-guardian django-pdb djangorestframework

#RUN pip install -r requirements.txt

#RUN conda develop /root/a2g2  # adds the root project dir to the pythonpath for this conda env
#RUN conda develop /root/a2g2/djangochem # adds the djangochem directory.  helpful for notebooks


# create a default bashrc
#RUN echo 'export PATH=/root/conda/bin:$PATH' >> /root/.bashrc

EXPOSE 9999
WORKDIR /a2g2/djangochem

#ENTRYPOINT /bin/bash --login
COPY /entrypoint.sh /
ENTRYPOINT ["/entrypoint.sh"]

#/bin/bash -c "pip install -e /a2g2/djangochem" &&
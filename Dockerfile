FROM osgeo/gdal:ubuntu-small-latest
RUN mkdir /pyfint
COPY . /pyfint
WORKDIR /pyfint
RUN apt-get update
RUN apt-get install -y python3-pip
RUN pip install -r requirements.txt
CMD python example_main.py

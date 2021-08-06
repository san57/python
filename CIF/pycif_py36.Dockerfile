# Build a Docker image that contains pyCIF and all requirements.
FROM python:3.6

RUN apt-get update \
	&& apt-get install -y libgdal-dev gdal-bin \
	&& apt-get autoclean && apt-get clean;

COPY . /tmp/CIF

RUN cd /tmp/CIF && pip install ".[dev, test]";

ENTRYPOINT ["pycif"]

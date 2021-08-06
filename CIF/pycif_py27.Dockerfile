# Build a Docker image that contains pyCIF and all requirements.
FROM python:2.7

RUN apt-get update \
	&& apt-get install -y libgdal-dev gdal-bin libeccodes-dev \
	&& apt-get autoclean && apt-get clean && pip install "numpy<1.17";

COPY . /tmp/CIF

RUN cd /tmp/CIF && pip install ".[dev, test]";

ENTRYPOINT ["pycif"]

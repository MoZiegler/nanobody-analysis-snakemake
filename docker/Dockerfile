FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
	gromacs \
	gromacs-data \
	python3 \
	python3-pip \
	python3-setuptools \
	python3-wheel \
	ca-certificates \
	&& gmx --version \
	&& rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN python3 -m pip install --no-cache-dir --upgrade pip \
	&& python3 -m pip install --no-cache-dir snakemake \
	&& python3 -m pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["snakemake", "--cores", "1"]

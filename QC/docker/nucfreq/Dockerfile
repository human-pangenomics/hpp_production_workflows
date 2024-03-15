FROM python:3.11-slim-buster

# Update package lists
RUN apt-get update && \
	apt-get install -y \
	wget \
	git && \
	apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/nucfreq

# Clone NucFreq repository
RUN git clone https://github.com/mrvollger/NucFreq.git .

# Install Python dependencies
RUN pip install --upgrade pip
RUN pip3 install numpy pysam pandas matplotlib seaborn

WORKDIR /data
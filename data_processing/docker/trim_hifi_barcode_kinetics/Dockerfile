FROM quay.io/biocontainers/pysam:0.23.3--py311hb456a96_0

MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL software="trim_hifi_barcode_kinetics" \ 
    container="trim_hifi_barcode_kinetics" \ 
    about.summary="Trim PacBio kinetics arrays (fi/fp/ri/rp) to match adapter‑clipped HiFi reads. This script may be necessary to run for files which were demultiplexed using lima <2.5. These files retained kinetics values for the barcodes but the sequence values did not resulting in a mismatch in their array lengths (and importantly a shift in the correct) position for the kinetics value in relation to the base prediction."

WORKDIR /opt
COPY scripts/trim_hifi_barcode_kinetics.py /opt/trim_hifi_barcode_kinetics.py

RUN chmod +x /opt/trim_hifi_barcode_kinetics.py


WORKDIR /data
FROM dfam/tetools:1.85
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL   \
    software="tetools" \ 
    container="tetools" \ 
    about.summary="dfam's tetools container with a rmToTrackHub.pl version that allows for avoiding bedtobigbed conversion."

## We want to convert *out and *align files to rmsk beds without creating
## bigbeds or an entire hub. This requires an old version of rmToTrackHub.pl
RUN curl -o \
    /opt/RepeatMasker/util/rmToTrackHub.pl \
    https://raw.githubusercontent.com/Dfam-consortium/RepeatMasker/fbc928138d652afdcafce6b73ffb7f7c84eccb9b/util/rmToTrackHub.pl
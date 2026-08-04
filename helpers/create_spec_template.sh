#!/usr/bin/env bash

cat << 'EOF' > spec_template.json
{
  "Software_version": "1.0.0",
  "subj": "",
  "session": "",
  "eddy_method": "default",
  "nstreamlines": "20M",
  "bidsdir": "<path/to/bidsfolder>",
  "outputdir": "<path/to/outputderivativesfolder>",
  "workdir": "<path/to/scratchworkdirectory>",
  "freesurferdir": "<path/to/existing/freesurferoutputdirectory",
  "nthreads": 8
}
EOF

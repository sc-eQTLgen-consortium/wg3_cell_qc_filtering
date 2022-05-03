#!/bin/bash
echo '# we need some more memory to do 
options(future.globals.maxSize = 120 * 1000 * 1024^2)
' > map_hierscpred.R
cat /map_hierscpred.R >> map_hierscpred.R
Rscript ./map_hierscpred.R "$@"


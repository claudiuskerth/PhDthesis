#!/bin/bash

for i in {1..100};
do
	realSFS -P 12 SAFs/ERY/ERY.FOLDED.saf.idx 2> /dev/null >> SFS/ERY/ERY.FOLDED.sfs.ml
	realSFS -P 12 SAFs/PAR/PAR.FOLDED.saf.idx 2> /dev/null >> SFS/PAR/PAR.FOLDED.sfs.ml
done

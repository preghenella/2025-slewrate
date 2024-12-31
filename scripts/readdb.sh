#! /usr/bin/env bash

if [ $# -ne 1 ]; then
    echo $0 [database]
    exit 1
fi
database=$(realpath $1)
localdir=$(dirname ${database})
echo ${localdir}

while read -r channel bcrconfig opmode deltathr deltathr2 vbias directory notes; do
    [[ $channel == \#* ]] && continue
    echo ${channel} ${bcrconfig} ${opmode} ${deltathr} ${deltathr2} ${vbias} ${directory} ${notes}

    subrun=$(basename ${directory})
    subrundir=${localdir}/subruns/${subrun}
    echo "local subrun directory: ${subrundir}"
    
done < ${database}

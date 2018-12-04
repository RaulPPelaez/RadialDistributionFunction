
nvcc --help |
    grep '\-\-gpu-code' -A1000 |
    grep -Po 'compute_\K[0-9]+' |
    sort |
    uniq |
    awk '{print "-gencode arch=compute_"$1",code=sm_"$1}' |
    paste -sd" " | tr '\n' ' '



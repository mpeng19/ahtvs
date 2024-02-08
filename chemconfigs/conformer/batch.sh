#!/bin/bash
#SBATCH --partition={{jobspec.details.partitions|join(",")}}
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 36:00:00
#SBATCH --mem-per-cpu=4096

{% for subdir in sub_dir_list %}
pushd {{subdir}}
source {{job_filename}}
popd
{% endfor %}
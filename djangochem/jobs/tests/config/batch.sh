{% for subdir in sub_dir_list %}
pushd {{subdir}}
source ./{{job_filename}}
popd
{% endfor %}
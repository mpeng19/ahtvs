#!/usr/bin/env nextflow

// params.lib_smi_uri="/home/denn/home2/ml/data/sf/Reaxys4_fps_uv.pickle"
params.lib_smi_uri="/home/denn/home2/ml/data/gp/reaxys4_s1_exp_theo_dftb.pickle"
lib_smi = file(params.lib_smi_uri)
smiles=file("sel.smi")
params.split_lines_smi=4

params.project = "demo"
params.tag = "d1"
params.dj_settings = "djangochem.settings.denn_tests"
params.chemconf_dir = "/home/denn/home/a2g2/chemconfigs"
params.conda_env = 'a2'
// process start_redis {

//     """
//     ~/bin/redis/redis-server /home/denn/home2/ml/data/redis-zinc1/zinc1.conf --dbfilename demo.db --unixsocket /tmp/redis_demo.sock
//     """

// }

process reset_db {
    
    output:
    val 1 into dbreset
    """
    . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
    export CPATH=""
    conda activate $params.conda_env
    python /home/denn/home/a2g2/djangochem/manage.py flush --noinput --settings=$params.dj_settings
    python /home/denn/home/a2g2/djangochem/manage.py migrate --noinput --settings=$params.dj_settings
    python /home/denn/home/a2g2/djangochem/manage.py scanconfigs $params.chemconf_dir --update --settings=$params.dj_settings
    """
}

// process random_smi {
//     input:
//     file lib_smi

//     output:
//     file "sel.smi" into smiles

//     """
//     . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
//     export CPATH=""
//     conda activate ml
//     python -m genchem.utils.smi_proc $lib_smi sel.smi 20
//     """

// }

process split_smi {
    
    input:
    file smiles

    output:
    file "reax4_*" into split_smiles mode flatten

    """
    split -l $params.split_lines_smi $smiles reax4_
    """
}


process add_smi {
    
    input:
    file smi from split_smiles
    val x from dbreset

    output:
    val 1 into (start_build, start_parse)
    """
    . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
    export CPATH=""
    conda activate $params.conda_env
    python /home/denn/home/a2g2/djangochem/manage.py addsmiles $params.project $smi -t $params.tag --settings=$params.dj_settings
    
    """
}

params.build_dir="/home/denn/home2/nf/act1//build"
params.tmp_dir="./tmp"
build_n=2

process build_jobs{
    input:
    val x from start_build.take(1)
    each y from 1..build_n

    output:
    stdout bld_job_out
    
    """
    . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
    export CPATH="" 
    conda activate $params.conda_env
    export DJANGO_SETTINGS_MODULE=$params.dj_settings
    export PIPE_TAG=$params.tag
    export PIPE_PROJECT=$params.project
    export PIPE_BUILD_DIR=$params.build_dir
    export PIPE_TMP_DIR=$params.tmp_dir
    export PIPE_BUILD_ITER=1
    python -m nextflow.nf_script_build_test 2>&1
    """
}

params.completed_dir="/home/denn/home2/nf/act1/complete"
params.arch_dir="/home/denn/home2/nf/act1/arch"

process parse_jobs{
    
    input:
    val x from start_parse.take(1)
    maxForks 3

    output:
    stdout parse_job_out
    
    """
    . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
    export CPATH="" 
    conda activate $params.conda_env
    export DJANGO_SETTINGS_MODULE=$params.dj_settings
    export PIPE_TAG=$params.tag
    export PIPE_PROJECT=$params.project
    export PIPE_COMPL_DIR=$params.completed_dir
    export PIPE_TMP_DIR=$params.tmp_dir
    export PIPE_ARCH_DIR=$params.arch_dir
    export PIPE_BUILD_ITER=1
    python -m nextflow.nf_script_parse
    """
}


bld_job_out.subscribe { println it }
parse_job_out.subscribe { println it }
// process calc_smi {
//     input:
//     file smiles

//     output:
//     file "calc_smi.csv" into calc_smiles

//     """
//     . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
//     export CPATH=""
//     conda activate ml
//     python -m genchem.utils.smi_calc $smiles calc_smi.csv
//     """

// }



// calc_smiles.subscribe {
//     println it
// }
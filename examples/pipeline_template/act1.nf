#!/usr/bin/env nextflow

// params.lib_smi_uri="/home/denn/home2/ml/data/sf/Reaxys4_fps_uv.pickle"
params.lib_smi_uri="/home/denn/home2/ml/data/gp/reaxys4_s1_exp_theo_dftb.pickle"
lib_smi = file(params.lib_smi_uri)

// process start_redis {

//     """
//     ~/bin/redis/redis-server /home/denn/home2/ml/data/redis-zinc1/zinc1.conf --dbfilename demo.db --unixsocket /tmp/redis_demo.sock
//     """

// }



process random_smi {
    input:
    file lib_smi

    output:
    file "sel.smi" into smiles

    """
    . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
    export CPATH=""
    conda activate ml
    python -m genchem.utils.smi_proc $lib_smi sel.smi
    """

}

process calc_smi {
    input:
    file smiles

    output:
    file "calc_smi.csv" into calc_smiles

    """
    . /home/denn/bin/miniconda3/etc/profile.d/conda.sh
    export CPATH=""
    conda activate ml
    python -m genchem.utils.smi_calc $smiles calc_smi.csv
    """

}


calc_smiles.subscribe {
    println it
}
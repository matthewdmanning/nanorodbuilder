#!/usr/bin/env bash
glob_pattern="$1"
mol2_dir="/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/"
antechamber_dir=${mol2_dir}/Renamed
mol2_backup_dir=Backup
overwrite=1
cd $mol2_dir
pwd
for mol2_file in *"$glob_pattern"*.mol2; do
    if [[ -f ${mol2_file} ]]; then
        antechamber_output_name=${mol2_file/'mol2'/'ac'}
        frcmod_file_name=${mol2_file/'mol2'/'frcmod'}
        if [[ ! -f ${antechamber_dir}/${antechamber_output_name} ]] || [[ ${overwrite} == 1 ]]; then
            echo "Creating ${antechamber_output_name}"
            atomtype -i ${mol2_file} -f mol2 -o ${antechamber_dir}/${antechamber_output_name} -p gaff2
            #cp ${mol2_file} ${mol2_backup_dir}
            #mv ${antechamber_dir}/${antechamber_output_name} ${mol2_file}
        fi
        if [[ ! -f frcmod/${frcmod_file_name} ]] || [[ ${overwrite} == 1 ]]; then
            echo "Creating ${frcmod_file_name}"
            parmchk2 -i ${mol2_file} -f mol2 -o frcmod/${frcmod_file_name} -s gaff2
            grep penalty frcmod/${frcmod_file_name}
            grep ATTN frcmod/${frcmod_file_name}
            #cp frcmod/${frcmod_file_name} ./${frcmod_file_name/'.ac'}.mol2
        fi
    fi
done

#/home/mdmannin/amber16/dat/leap/parm/AuNP.frcmod
#/bin/bash
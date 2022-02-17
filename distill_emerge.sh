
# python ./make_emerge_modules.py
# python /home/rman/projects/emerge/DRAM_tools/DRAM_tools/distill_genomes.py --help

python3 ./emerge_distill_tools/distill_genomes1_2.py \
        -i ./dram_data/all_emerge_annotations.tsv \
        -o emerge_distilled_output \
        --replace_forms \
        --custom_distillate ./custom_input_modules/EMERGE_distillate_module.tsv \
        --custom_pathways ./custom_input_modules/EMERGE_pathways_module.tsv 


#DONE

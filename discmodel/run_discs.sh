bgadd -L 200 /discs
for i in $(seq 0 199)
do bsub -g /discs -o console_output/beastLV_out_"$i".txt -e console_output/beastLV_err_"$i".txt python discmodel/discs.py -jobi "$i" -N 1 --re_run
done

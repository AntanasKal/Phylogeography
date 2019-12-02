bgadd -L 200 /LV_phyrex 
for i in $(seq 0 199)
do bsub -g /LV_phyrex -J LV"$i" -o console_output/phyrex_out_"$i".txt -e console_output/phyrex_err_"$i".txt python launch_phyrex.py -index $i -sample_index 4
done

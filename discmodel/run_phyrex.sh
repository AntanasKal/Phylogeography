bgadd -L 200 /LV_phyrex 
for i in $(seq 106 106)
do bsub -g /LV_phyrex -o console_output/phyrex_out_"$i".txt -e console_output/phyrex_err_"$i".txt ./phyrex --xml=output/phyrex/LV/phyrex_input/phyrex"$i".xml
done

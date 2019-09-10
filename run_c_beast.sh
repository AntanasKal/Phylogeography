for i in $(seq 0 99)
do bgadd -L 10 /c_beast/sim"$i"
bsub -g /c_beast/sim"$i" -o console_output/phyrex_out1_"$i".txt -e console_output/phyrex_err1_"$i".txt python launch_corrected_beast.py -index $i -sample_index 1
bsub -g /c_beast/sim"$i" -o console_output/phyrex_out2_"$i".txt -e console_output/phyrex_err2_"$i".txt python launch_corrected_beast.py -index $i -sample_index 2
bsub -g /c_beast/sim"$i" -o console_output/phyrex_out3_"$i".txt -e console_output/phyrex_err3_"$i".txt python launch_corrected_beast.py -index $i -sample_index 3
bsub -g /c_beast/sim"$i" -o console_output/phyrex_out4_"$i".txt -e console_output/phyrex_err4_"$i".txt python launch_corrected_beast.py -index $i -sample_index 4
done

for i in $(seq 0 99)
do bgadd -L 10 /c_beast/sim"$i"
bsub -g /c_beast/sim"$i" -o console_output/c_beast_out1_"$i".txt -e console_output/c_beast_err1_"$i".txt python launch_corrected_beast.py -index $i -sample_index 1 --re_run
bsub -g /c_beast/sim"$i" -o console_output/c_beast_out2_"$i".txt -e console_output/c_beast_err2_"$i".txt python launch_corrected_beast.py -index $i -sample_index 2 --re_run
bsub -g /c_beast/sim"$i" -o console_output/c_beast_out3_"$i".txt -e console_output/c_beast_err3_"$i".txt python launch_corrected_beast.py -index $i -sample_index 3 --re_run
bsub -g /c_beast/sim"$i" -o console_output/c_beast_out4_"$i".txt -e console_output/c_beast_err4_"$i".txt python launch_corrected_beast.py -index $i -sample_index 4 --re_run
done

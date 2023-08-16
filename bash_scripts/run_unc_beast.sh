bgadd -L 400 /unc_beast
for i in $(seq 0 99)
do bsub -g /unc_beast -o console_output/unc_beast_out1_"$i".txt -e console_output/unc_beast_err1_"$i".txt python launch_uncorrected_beast.py -index $i -sample_index 1 --re_run
bsub -g /unc_beast -o console_output/unc_beast_out2_"$i".txt -e console_output/unc_beast_err2_"$i".txt python launch_uncorrected_beast.py -index $i -sample_index 2 --re_run
bsub -g /unc_beast -o console_output/unc_beast_out3_"$i".txt -e console_output/unc_beast_err3_"$i".txt python launch_uncorrected_beast.py -index $i -sample_index 3 --re_run
bsub -g /unc_beast -o console_output/unc_beast_out4_"$i".txt -e console_output/unc_beast_err4_"$i".txt python launch_uncorrected_beast.py -index $i -sample_index 4 --re_run
done

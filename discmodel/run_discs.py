bgadd -L 200 /discs
for i in $(seq 107 107)
do bsub -g /discs -o console_output2/out_"$i".txt -e console_output2/err_"$i".txt python disc2.py -jobi "$i" -N 1
done

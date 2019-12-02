for i in $(seq 1 1000)
do bsub -o console_output/out_"$i".txt -e console_output/err_"$i".txt python disc_simulation.py -jobi "$i" -N 10
done

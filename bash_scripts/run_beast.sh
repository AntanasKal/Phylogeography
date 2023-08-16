for i in $(seq 0 99)
do bsub -o console_output/out_"$i".txt -e console_output/err_"$i".txt python simulation.py -jobi "$i" -N 1 -dims 2 -treetype yule -ntips 1000 --c_beast --re_run --linux -mcmc 10000
done

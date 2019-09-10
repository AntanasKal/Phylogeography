for i in $(seq 0 99)
do bsub -o console_output/out_"$i".txt -e console_output/err_"$i".txt python simulation.py -jobi "$i" -N 1 -dims 2 -treetype nuc -ntipspp 100 -nps 20 --linux -mcmc 10000
done

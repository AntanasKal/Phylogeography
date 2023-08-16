for j in $(seq 1 5)
do for i in $(seq 0 99)
do bsub -o console_output/out_extra_"$j"_"$i".txt -e console_output/err_extra_"$j"_"$i".txt python simulation_extra.py -jobi "$i" -N 1 -dims 2 -mcmc 10000 -biasIntensity -gradient "$j"
done
done

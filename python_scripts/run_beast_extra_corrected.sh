for j in $(seq 1 5)
do for i in $(seq 0 99)
do bsub -o console_output/out_extraC_"$j"_"$i".txt -e console_output/err_extraC_"$j"_"$i".txt python simulation_extra.py -jobi "$i" -N 1 -dims 2 -extraSamples -gradient "$j"
done
done

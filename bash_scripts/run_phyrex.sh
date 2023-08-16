bgadd -L 400 /yule/phyrex
for i in $(seq 0 99)
do bsub -g /yule/phyrex -o console_output/phyrex_out1_"$i".txt -e console_output/phyrex_err1_"$i".txt /nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/sampled1/phyrex_input/phyrex"$i".xml
bsub -g /yule/phyrex -o console_output/phyrex_out2_"$i".txt -e console_output/phyrex_err2_"$i".txt /nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/sampled2/phyrex_input/phyrex"$i".xml
bsub -g /yule/phyrex -o console_output/phyrex_out3_"$i".txt -e console_output/phyrex_err3_"$i".txt /nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/sampled3/phyrex_input/phyrex"$i".xml
bsub -g /yule/phyrex -o console_output/phyrex_out4_"$i".txt -e console_output/phyrex_err4_"$i".txt /nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/sampled4/phyrex_input/phyrex"$i".xml
done


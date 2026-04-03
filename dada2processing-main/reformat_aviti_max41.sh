for i in *.fastq.gz; do echo "/home/umii/goul0109/bbmap/reformat.sh in=$i out=../reformat/${i//.fastq.gz/.reformat.fastq.gz} mincalledquality=2 maxcalledquality=41 qin=33" >> run_bbreformat.cmd; done
chmod +x run_bbreformat.cmd
./run_bbreformat.cmd
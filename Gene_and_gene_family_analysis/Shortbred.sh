{\rtf1\ansi\ansicpg936\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 ArialMT;}
{\colortbl;\red255\green255\blue255;\red26\green26\blue26;\red255\green255\blue255;\red16\green60\blue192;
}
{\*\expandedcolortbl;;\cssrgb\c13333\c13333\c13333;\cssrgb\c100000\c100000\c100000;\cssrgb\c6667\c33333\c80000;
}
\paperw11900\paperh16840\margl1440\margr1440\vieww38200\viewh21040\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs36 \cf2 \cb3 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 ###############################\
Step1: identify\cb1 \
\cb3 #!/bin/bash\cb1 \
\cb3 #SBATCH --job-name=shortbred\cb1 \
\cb3 #SBATCH --error=shortbred.err\cb1 \
\cb3 #SBATCH --output=shortbred.out\cb1 \
\cb3 #SBATCH --mem=60gb\cb1 \
\cb3 #SBATCH --time=72:00:00\cb1 \
\cb3 #SBATCH --cpus-per-task=8\cb1 \
\cb3 ml Anaconda3/2022.05\cb1 \
\cb3 conda activate /scratch/hb-tifn/condas/conda_assemblers\cb1 \
\
\cb3 shortbred_identify.py --goi /scratch/hb-fu/yzhang/Progenome1_bakta_flanking1_MAGgenes.faa --ref /scratch/p303998/SV_MWAS/Traitar/uniref90.fasta --markers Two_methods_109SVs_1375genes.faa --tmp tmp_markers --usearch /scratch/hb-tifn/tools/GMH_pipeline/utils/usearch10.0.240_i86linux32 --threads 8 --blastp /scratch/hb-tifn/condas/conda_assemblers/bin/blastp --cdhit /scratch/hb-tifn/condas/conda_assemblers/bin/cd-hit --makeblastdb /scratch/hb-tifn/condas/conda_assemblers/bin/makeblastdb\cb1 \
\
\cb3 Step2: Quantify\cb1 \
\
\cb3 ##### step_1: extract samples' names\cb1 \
\cb3 for sample in *_1.fq.gz\cb1 \
\cb3 do i=$(echo $\{sample%_1.fq.gz\})\cb1 \
\
\cb3 ##### step_2: basic settings\cb1 \
\cb3 echo "#!/bin/bash" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "#SBATCH --job-name=$i" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "#SBATCH --error=$i.shortbred.err" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "#SBATCH --output=$i.shortbred.out" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "#SBATCH --mem=60gb" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "#SBATCH --time=20:00:00" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "#SBATCH --cpus-per-task=4" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\
\cb3 echo "zcat $i\\_1.fq.gz > $i.merge.fastq; zcat $i\\_2.fq.gz >> $i.merge.fastq" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\
\cb3 echo "ml Anaconda3/2022.05" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "conda activate /scratch/hb-tifn/condas/conda_assemblers" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\
\cb3 echo "shortbred_quantify.py --markers /scratch/hb-fu/yzhang/Two_methods_109SVs_1375genes.faa --wgs $i.merge.fastq --threads 16 --results $i.noclean.txt --tmp $i/tmp --usearch /scratch/hb-tifn/tools/GMH_pipeline/utils/usearch10.0.240_i86linux32" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "echo\'a0shortbred\'a0finished" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "rm $i.merge.fastq" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\cb3 echo "rm $i/ -rf" >> ./${\field{\*\fldinst{HYPERLINK "http://i.shortbred.sh/"}}{\fldrslt \cf4 \ul \ulc4 \strokec4 i.shortbred.sh}}\cb1 \
\
\cb3 done\cb1 \
\cb3 #$i\\_1.fq.gz $i\\_2.fq.gz\
}
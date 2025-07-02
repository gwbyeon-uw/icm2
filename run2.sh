for file in `ls fastq_use/*R1*.fastq`; 
do 
	name=${file%.fastq}; 
	name=$(basename $name | cut -d"_" -f2); 
	r2=$(echo $file | sed s/_R1_/_R2_/g); 
	echo $name $file $r2;
	
	source E06/E06_"$name".config

	mkdir -p out/$name
	if [ "$SPLIT" == "0" ]; then
		mkdir -p out/$name/split0
		awk ' $0~/^>/ { header=$0; getline; print header; print "'$CONST5_0'"$0"'$CONST3_0'" } ' $REF_0 > out/$name/split0/$name"_split0.fa" #To build FASTA with constant regions from config added

		./p0.sh $file $r2 out/$name $name false $MERGE_0
		./p1.sh out/$name/$name"_split0.fastq" out/$name/split0 $name"_split0" true $UMI_0 out/$name/split0/$name"_split0.fa" $CONST5_0 $CONST3_0 ${#CONST5_0} ${#CONST3_0}
	else
		mkdir -p out/$name/split1
		mkdir -p out/$name/split2

		awk ' $0~/^>/ { header=$0; getline; print header; print "'$CONST5_1'"$0"'$CONST3_1'" } ' $REF_1 > out/$name/split1/$name"_split1.fa" #To build FASTA with constant regions from config added
		awk ' $0~/^>/ { header=$0; getline; print header; print "'$CONST5_2'"$0"'$CONST3_2'" } ' $REF_2 > out/$name/split2/$name"_split2.fa"

	        ./p0.sh $file $r2 out/$name $name true $MERGE_1 
		./p1.sh out/$name/$name"_split1.fastq" out/$name/split1 $name"_split1" true $UMI_1 out/$name/split1/$name"_split1.fa" $CONST5_1 $CONST3_1 ${#CONST5_1} ${#CONST3_1}
        	./p1.sh out/$name/$name"_split2.fastq" out/$name/split2 $name"_split2" true $UMI_2 out/$name/split2/$name"_split2.fa" $CONST5_2 $CONST3_2 ${#CONST5_2} ${#CONST3_2}
	fi
done

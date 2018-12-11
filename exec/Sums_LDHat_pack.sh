#! /bin/sh
seqName=$1
nn=$2
lpart=$3
pathLDhat=$4
jobId=$5
out=$6
curr=$(pwd)
# sed -i "1i $nn $lpart 1" $seqName

awk -v n=1 -v s="$nn $lpart 1" 'NR == n {print s} {print}' $seqName > ${curr}/"FastaTemp_${out}${jobId}.fa"
mv ${curr}/"FastaTemp_${out}${jobId}.fa" $seqName
mkdir -p "${out}${jobId}" && cd "${out}${jobId}"
$pathLDhat/convert -seq $seqName | tail -n 7 >"temp_${out}${jobId}.txt"
cat "temp_${out}${jobId}.txt" >>../"Sums_part_main_${out}${jobId}.txt"

# cd $pathLDhat/lk_files
# if [ ! -f "lk_n${nn}_t${th}_rh100_npts201.txt" ]; then
# echo "Lookup table is calculated"
# $pathLDhat/complete -n $nn -rhomax 100 -n_pts 201 -theta $th -split $cor
# mv new_lk.txt "lk_n${nn}_t${th}_rh100_npts201.txt"
# fi


# pathLDhat=$(locate LDhat-master | head -n 1)
# pathLooks=$(locate "lk_n${nn}_t${th}_rh100_npts201.txt" | head -n 1)
#### make here subdirectory for each part and then write files somewhere outside (resLDH...) ####
#cd $curr/${out}${jobId}
#echo 0 | $pathLDhat/pairwise -seq sites.txt -loc locs.txt -lk $pathLDhat/lk_files/"lk_n${nn}_t${th}_rh100_npts201.txt" -prefix ${out}${jobId} >"test${out}.txt"
# [ -f outfile.txt ] && echo "File exist" || echo "File does not exist"
#head -n 5 "${out}${jobId}outfile.txt" | tail -n 1 >"${out}${jobId}temp.txt"
#cat "${out}${jobId}temp.txt" >>../"resLDHats_pairwise_main_${out}${jobId}.txt"

cd .. && rm -rf ${out}${jobId}

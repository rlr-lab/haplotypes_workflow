#! /bin/bash
set -eu
### required arguments
file_sam=""
file_ref=""
### optional arguments
file_path='./result'
prefix="rvhaplo"
mq=0
thread=8
error_rate=0.1
signi_level=0.05
cond_pro=0.65
fre_snv=0.8
num_read_1=10
num_read_2=5
gap=15
smallest_snv=20
only_snv=0
ovlap_read=5
weight_read=0.85
mcl_inflation=2
lar_cluster=50
ovlap_cluster=10
depth=5
weight_cluster=0.8
abundance=0.005
s_pos=1
e_pos=10000000000
sub_graph=1


function help_info() {
	echo "Usage: ${0-} -i alignment.sam -r ref_genome.fasta [options]"
	echo ""
	echo "RVHaplo: Reconstructing viral haplotypes using long reads"
	echo ""
	echo "Author: Dehan CAI"
	echo "Date:   May 2022"
	echo "Version 2: Support mutli-thread processing; Use a C package of MCL; Cost less memory   "
	echo ""
	echo "    -i | --input:                     alignment file (sam)"
	echo "    -r | --refernece:                 reference genome (fasta)"
	echo ""
	echo "    Options:"
	echo "    -o  | --out:                      Path where to output the results. (default:./result)"
	echo "    -p  | --prefix STR:               Prefix of output file. (default: rvhaplo)"
	echo "    -t  | --thread INT:               Number of CPU cores for multiprocessing. (default:8)"
	echo "    -e  | --error_rate FLOAT:         Sequencing error rate. (default: 0.1)"
	echo "    -mq | --map_qual INT:             Smallest mapping quality for reads . (default:0)"
	echo "    -s  | --signi_level FLOAT:        Significance level for binomial tests. (default: 0.05)"
	echo "    -c  | --cond_pro FLOAT:           A threshold of the maximum conditional probability for SNV sites. (default: 0.65)"
	echo "    -f  | --fre_snv FLOAT:            The most dominant base' frequency at a to-be-verified site should >= fre_snv. (default: 0.80)"
	echo "    -n1 | --num_read_1 INT:           Minimum # of reads for calculating the conditional probability given one conditional site. (default:10)"
	echo "    -n2 | --num_read_2 INT:           Minimum # of reads for calculating the conditional probability given more than one conditional sites. (default: 5)"
	echo "    -g  | --gap INT:                  Minimum length of gap between SNV sites for calculating the conditional probability. (default:15)"
	echo "    -ss | --smallest_snv INT:         Minimum # of SNV sites for haplotype construction. (default:20)"
	echo "    -os | --only_snv (0 or 1) :       Only output the SNV sites without running the haplotype reconstruction part. (default: 0)"
	echo "    -or | --overlap_read INT:         Minimum length of overlap for creating edges between two read in the read graph. (default: 5)"
	echo "    -wr | --weight_read FLOAT:        Minimum weights of edges in the read graph. (default:0.85)"
	echo "    -sg | --sub_graph INT:            Number of subgraphs to run MCL (default:1)"
	echo "    -m  | --mcl_inflation FLOAT:      Inflation of MCL algorithm. (default:2)"
	echo "    -l  | --lar_cluster INT:          A threshold for seperating clusters into two groups based on sizes of clusters. (default:50)"
	echo "    -oc | --overlap_cluster INT:      A parameter related to the minimum overlap between consensus sequences. (default:10) "
	echo "    -d  | --depth INT:                Depth limitation for consensus sequences generated from clusters. (default:5) "
	echo "    -wc | --weight_cluster FLOAT:     Minimum weights between clusters in the hierarchical clustering. (default: 0.8)"
	echo "    -sp | --start_pos INT:            Starting position for generating consensus sequences (default: 1)"
	echo "    -ep | --end_pos INT:              Ending position for generating consensus sequences. (default: 1e10)"
	echo "    -a  | --abundance FLOAT:          A threshold for filtering low-abundance haplotypes. (default: 0.005)"
	echo "    -h  | --help :                    Print help message."
	echo ""
	echo "    For further details of above arguments, please refer to https://github.com/dhcai21/RVHaplo   "
	echo ""
	exit 0
}

if [[ "${1-}" == "" ]];then
	help_info
	exit 1
fi

while [[ "${1-}" != "" ]]; do
	case "${1-}" in
		-h | --help ) ## print help message
		help_info
		exit 1
		;;
		-i | --input ) ### input sam file
		case "${2-}" in 
		"" )
			echo "Error: no sam file as input"
			exit 1
			;;
		*)
			file_sam="${2-}"
			if [[ "${file_sam:0:1}" == "-" ]]
			then
				echo "Error: no sam file as input"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-r | --ref_genome) ### input reference genome
		case "${2-}" in 
		"")
			echo "Error: no fasta file as input"
			exit 1
			;;
		*)
			file_ref="${2-}"
			if [[ ""${file_ref:0:1}"" == "-" ]]
			then
				echo "Error: no fasta file as input"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-o | --out )  ### output path
		case "${2-}" in 
		"" )
			echo "Error: no output path"
			exit 1
			;;
		*)
			file_path="${2-}"
			if [[ "${file_sam:0:1}" == "-" ]]
			then
				echo "Error: no output path"
				exit 1
			fi
			shift 2
			;;
		esac
		;;
		-p | --prefix )  ### prefix
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			prefix="${2-}"
			shift 2
			;;
		esac
		;;
		-mq | --map_qual )  ### mapping quality
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			mq="${2-}"
			shift 2
			;;
		esac
		;;
		-t | --thread )  ### threads
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			thread="${2-}"
			shift 2
			;;
		esac
		;;
		-e | --error_rate )  ### error_rate
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			error_rate="${2-}"
			shift 2
			;;
		esac
		;;
		-s | --signi_level )  ### significance_level
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			signi_level="${2-}"
			shift 2
			;;
		esac
		;;
		-c | --cond_pro )  ### conditional_probability
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			cond_pro="${2-}"
			shift 2
			;;
		esac
		;;
		-f | --fre_snv )  ### determine the set of to-be-verified SNV sites
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			fre_snv="${2-}"
			shift 2
			;;
		esac
		;;
		-n1 | --num_read_1 )  ### number of reads for p(ai|aj)
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			num_read_1="${2-}"
			shift 2
			;;
		esac
		;;
		-n2 | --num_read_2 )  ### number of reads for p(ai|aj1,aj2,...,ajp)
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			num_read_2="${2-}"
			shift 2
			;;
		esac
		;;
		-g | --gap )  ### Minimum distance between SNV sites
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			gap="${2-}"
			shift 2
			;;
		esac
		;;
		-ss | --smallest_snv )  ### Minimum number of SNV sites for haplotype reconstruction
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			smallest_snv="${2-}"
			shift 2
			;;
		esac
		;;
		-os | --only_snv )  ### Only output the SNV sites without running the haplotype reconstruction part.
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			only_snv="${2-}"
			shift 2
			;;
		esac
		;;
		-or | --ovlap_read )  ### overlap_read
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			ovlap_read="${2-}"
			shift 2
			;;
		esac
		;;
		-wr | --weight_read )  ### weight_read
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			weight_read="${2-}"
			shift 2
			;;
		esac
		;;
		-sg | --sub_graph )  ### Number of subgraphs
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			sub_graph="${2-}"
			shift 2
			;;
		esac
		;;
		-m | --mcl_inflation )  ### inflation of MCL
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			mcl_inflation="${2-}"
			shift 2
			;;
		esac
		;;
		-oc | --ovlap_cluster )  ### overlap_cluster
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			ovlap_cluster="${2-}"
			shift 2
			;;
		esac
		;;
		-wc | --weight_cluster )  ### weight_cluster
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			weight_cluster="${2-}"
			shift 2
			;;
		esac
		;;
		-d | --depth )  ### depth limitation
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			depth="${2-}"
			shift 2
			;;
		esac
		;;
		-l | --lar_cluster )  ### large cluster size
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			lar_cluster="${2-}"
			shift 2
			;;
		esac
		;;
		-sp | --start_pos )  ### start_pos
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			s_pos="${2-}"
			shift 2
			;;
		esac
		;;
		-ep | --end_pos )  ### end_pos
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			e_pos="${2-}"
			shift 2
			;;
		esac
		;;
		-a | --abundance )  ### smallest abundance
		case "${2-}" in 
		"" )
			echo "Error: no input for ${1-}"
			exit 1
			;;
		*)
			abundance="${2-}"
			shift 2
			;;
		esac
		;;
		*)
			echo "Error: unknow parameter ${1-}"
			exit 1
	esac
done

if [[ "$file_sam" == "" ]];then
	echo "Error: no sam file input"
	exit 1
fi

if [[ "$file_ref" == "" ]];then
	echo "Error: no reference genome input"
	exit 1
fi

if [[ ${file_path:0-1:1} == "/" ]];then
	path_len=`expr ${#file_path}-1`
	file_prefix=$file_path$prefix
	file_path=${file_path:0:path_len}
else
	file_prefix=$file_path"/"$prefix
fi

current_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P)
workflow_path=$(dirname $current_path)
if [ -z "$workflow_path" ]; then
	workflow_path="/home/lzh8485/haplotypes_workflow"
fi
##########  count nucleotide occurrence  ##########
echo "count nucleotide occurrence"
rm -rf $file_path"/alignment"
mkdir -p $file_path"/alignment"
file_len=`expr ${#file_sam}-4`
unique_sam=$file_path"/alignment/"$prefix".sam"
samtools view -h -F 0x900 -q $mq $file_sam > $unique_sam
file_bam=$file_path"/alignment/"$prefix".bam"
samtools view -b -S $unique_sam > $file_bam
rm $unique_sam
file_bam_sorted=$file_path"/alignment/"$prefix"_sorted.bam"
samtools sort $file_bam -o $file_bam_sorted
samtools index $file_bam_sorted
file_acgt=$file_prefix"_acgt.txt"
python ${workflow_path}/RVHaplo/count_frequency.py $file_bam_sorted $file_acgt

########## two binomial tests  ##########
echo "SNV detection"
file_snv=$file_prefix"_snv.txt"
python ${workflow_path}/RVHaplo/two_binomial.py $error_rate $signi_level $file_acgt $file_snv $thread $s_pos $e_pos

## judge number of detected SNV sites
size="$(wc -l <"$file_snv")"
#size="${size:0-1:1}"
if (( $size != 0 ));then
	python ${workflow_path}/RVHaplo/out_haplotypes.py $file_prefix"_clusters.pickle" $file_bam_sorted $file_path $file_acgt 1 $file_prefix"_consensus.fasta" $s_pos $e_pos $workflow_path
	python ${workflow_path}/RVHaplo/extract_reads.py $file_path $prefix 1
	python ${workflow_path}/RVHaplo/run_medaka.py $file_path $prefix 1
	rm -rf $file_path"/medaka"
	exit 0
fi

## maximum conditional probability and construct reads graph
python ${workflow_path}/RVHaplo/mcp_read_graph.py $file_bam_sorted $file_snv $cond_pro $smallest_snv $num_read_1 $num_read_2 $gap \
	$weight_read $ovlap_read $file_prefix $fre_snv $thread $only_snv $sub_graph

## judge number of detected SNV sites
size="$(wc -l <"$file_snv")"
#size="${size:0-1:1}"
if (( $size != 0 ));then
	python ${workflow_path}/RVHaplo/out_haplotypes.py $file_prefix"_clusters.pickle" $file_bam_sorted $file_path $file_acgt 1 $file_prefix"_consensus.fasta" $s_pos $e_pos $workflow_path
	python ${workflow_path}/RVHaplo/extract_reads.py $file_path $prefix 1
	python ${workflow_path}/RVHaplo/run_medaka.py $file_path $prefix 1
	rm -rf $file_path"/medaka"
	exit 0
fi

if (( $only_snv != 0 ));then
	exit 0
fi

## check the number of reads with overlaps
if (( $sub_graph != 1 )); then
	size="$(wc -l <"$file_prefix"0_reads_graph.txt)"
else
	size="$(wc -l <"$file_prefix"_reads_graph.txt)"
fi
#size="${size:0-1:1}"
if (( $size == 0 ));then
	echo "Not enough reads with overlaps"
	exit 0
fi

# MCL clustering
echo "MCL clustering"
python ${workflow_path}/RVHaplo/run_mcl.py "$file_prefix" "$thread" "$mcl_inflation" "$sub_graph"

## hierarchical clustering
echo "hierarchical clustering"
python ${workflow_path}/RVHaplo/hierarchical_cluster.py $file_prefix"_matrix.pickle" $lar_cluster $depth \
	$ovlap_cluster $weight_cluster $abundance $file_prefix

## reconstruct haplotypes
rm -rf $file_path"/clusters"
mkdir -p $file_path"/clusters"

echo "haplotypes reconstruction"

python ${workflow_path}/RVHaplo/out_haplotypes.py $file_prefix"_clusters.pickle" $file_bam_sorted $file_path $file_acgt x $file_prefix"_consensus.fasta" $s_pos $e_pos $workflow_path

echo "haplotypes polishment (medaka)"
python ${workflow_path}/RVHaplo/extract_reads.py $file_path $prefix x
python ${workflow_path}/RVHaplo/run_medaka.py $file_path $prefix x

rm $file_prefix"_matrix.pickle"
rm $file_prefix"_reads_cluster.txt"
rm $file_prefix"_clusters.pickle"
rm -rf $file_path"/medaka"
echo "complete reconstructing haplotypes"

exit 0

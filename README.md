# Loci Networks from LD-clumping data
Finding number of "loci" using clumping data from LD-clumping

## Usage

./graph_network_analysis.py --types $clump1 $clump2 ... --files $clump1_file $clump2_file ... --analysis $clump1_$clump2_analysis --summary_dir $dir

## Output

There will be an output file with the name specified as "--analysis" flag in "--summary_dir" flag showing the whole VENN diagram of signals;

For the example given in usage part, the summary file will look like;

$clump1 : 10 -> means there are 10 exclusive signals found in clump1

$clump2  : 15 -> means there are 15 exclusive signals found in clump2

$clump1_$clump2: 100 -> means there are 100 signals both found in clump1 and clump2 

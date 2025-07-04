# GraphGerbil
Ploidy-aware genome phasing by graph partitioning, implemented in golang.

## To do list 🐁

- [ ] Refactor to have a parameter which will return the haplotype fastas as partitioned. Currently, I am only returning general stats of the partition process, not its actual outcome.
- [ ] Refactor to permit running using Minimisers only, or Kmers only - as desired by users. Current implementation can do Syncmers only or BUSCOs only.

## Premise

<p align="center">
  
<img width="749" alt="Screenshot 2025-07-04 at 15 02 03" src="https://github.com/user-attachments/assets/89a3be9f-9992-4a91-ad37-40e40efcfcc4" />

</p>

### How to run

Currently, two modes are supported. In syncmers mode, the tool will generate syncmers for the fasta file, and use those as markers to partition the genome. In busco mode, the tool will only parition the sequences which have BUSCOs on them.

```
Usage of ../GraphGerbil2:

  -Mode string
    	Phasing mode to use: 'syncmers' or 'busco'
  -Fasta string
    	(non-BUSCO modes only) Unphased fasta file
  -Fai string
    	Samtools index .fai file for the unphased fasta
  -K int
    	(non-BUSCO modes only) Kmer length, default 31
  -S int
    	(non-BUSCO modes only) Syncmer S, default 11
  -T int
    	(non-BUSCO modes only) Syncmer T, default 1
  -BuscoFullTable string
    	BUSCO full_table.tsv file for the unphased fasta
  -Max int
    	(non-BUSCO modes only) Maximum genome-wide count of markers used in phasing, default is Ploidy
  -Min int
    	(non-BUSCO modes only) Minimum genome-wide count of markers used in phasing, default is 2
  -Ploidy int
    	Ploidy, or number of bins you want to partition your genome into
  -R int
    	(non-BUSCO modes only) Number of markers to subsample when building the graph, default is to use all markers available (default -1)
  -cores int
    	Number of cores you want the programme to use, default is all available
  -DuplicationThresholds string
    	Starting thresholds to test, please provide as comma-separated list, e.g.: 1,2,3 (default "1,10,20,30,40,50")
```

### Some testing results


1. Carpinus betulus

Syncmer mode on dhCarBetu1 (Octoploid), with expected haploid genome size of ~400 Mb and ~30% BUSCO duplication. See: https://tolqc.cog.sanger.ac.uk/darwin/dicots/Carpinus_betulus/

```
Launching SyncmerMode on dhCarBetu1.fasta with 40 cores …
Time elapsed to process Samtools fai file: 3.311478ms
Time elapsed to process BUSCO files: 14.662284ms
4210 sequences in total, 1016 of those have BUSCOs on them 
Time elapsed to generate Syncmers for the whole genome: 45.397773931s
Time elapsed to generate SyncmerMap from SyncmerList: 2m3.294190401s
Time elapsed to filter syncmers: 3.871111993s
Time elapsed to generate SyncmerSequenceMap from Filtered SyncmerMap: 2.523574093s
3216 sequences in total are represented by 4190751 valid syncmers 
Time elapsed to sample 500000 syncmers: 542.710922ms
2288 sequences in total are represented by sampled syncmers 
Time elapsed to build graph: 2.344621931s
---------------------------------------------------
Threshold  Bin   Seqs  Span       C(%)   S(%)   D(%)  
---------------------------------------------------
10.0       1     1675  334000429  61.66  42.22  19.45 
10.0       2     242   428329344  82.52  51.04  31.47 
10.0       3     126   427761432  81.02  51.60  29.42 
10.0       4     87    450818800  81.41  49.72  31.68 
10.0       5     64    408883121  77.19  53.39  23.80 
10.0       6     45    379572885  71.56  52.88  18.68 
10.0       7     33    373018908  65.42  48.27  17.14 
10.0       8     12    165878364  36.50  34.54  1.96  
Unassigned       4     51626333   12.58  12.49  0.09
---------------------------------------------------
```

2. Allolobophora molleri

Syncmer mode on whApoMoll1 (tetraploid), with expected haploid genome size of ~1.5 Gb and ~15% BUSCO duplication. See: https://tolqc.cog.sanger.ac.uk/erga-bge/annelids/Allolobophora_molleri/

```

```


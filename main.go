package main

import (
	"GraphGerbil2/BUSCO"
	"GraphGerbil2/FASTK"
	"GraphGerbil2/Fai"
	"GraphGerbil2/Graph"
	"GraphGerbil2/Kmers"
	"GraphGerbil2/Sequence"
	"flag"
	"fmt"
	"os"
	"time"
)

func BUSCOMode(FaiPath string, BUSCOFullTable string, NumberOfBins int, ThresholdsToTest []float32) {

	// === STEP 1: Read BUSCO Full Table ===
	SequenceMap, SequenceHashMap, _, BUSCOMap := BUSCO.ReadBUSCOFullTable(BUSCOFullTable)
	TotalBUSCONum := len(BUSCOMap)

	// === STEP 2: Read sequence lengths ===
	SequenceLengthMap := Fai.ReadFaiFile(FaiPath, SequenceHashMap)

	// === STEP 3: Construct graph from BUSCOs ===
	BufferGraph, graph := Graph.GenerateGraph(BUSCOMap)
	SortedSequences := Graph.SortAllSequencesbyConnectedness(graph)

	fmt.Println("Threshold\tBin\tSeqs\tSpan\tS(%)\tD(%)")

	for _, Threshold := range ThresholdsToTest {
		// === STEP 4: Bin sequences using graph ===
		Bins, UnassignedSequences := Graph.PartitionGraph(
			BufferGraph,
			graph,
			SequenceMap,
			SortedSequences,
			NumberOfBins,
			Threshold,
			TotalBUSCONum,
			SequenceLengthMap,
		)

		// === STEP 5: Score and report each bin using BUSCOs ===
		for i, Bin := range Bins {
			fmt.Printf("%.1f\t%d\t%d\t%d\t%.2f\t%.2f\n",
				Threshold,
				i+1,
				Bin.NumofSequences,
				Bin.AssemblySpan,
				Bin.Score.ActualSingleCopyPercentage,
				Bin.Score.ActualDuplicatePercentage,
			)
		}

		// === STEP 6: Score unassigned sequences ===
		var UnassignedSequencesBUSCOs []uint64
		for _, Sequence := range UnassignedSequences {
			UnassignedSequencesBUSCOs = append(UnassignedSequencesBUSCOs, SequenceMap[Sequence].Markers...)
		}
		UnassignedSequencesScore := Sequence.ScoreMarkers(UnassignedSequencesBUSCOs, TotalBUSCONum)
		UnassignedSequencesSpan := Fai.MeasureAssemblySpan(UnassignedSequences, SequenceLengthMap)

		fmt.Printf("Unassigned\t%d\t%d\t%.2f\t%.2f\n",
			len(UnassignedSequences),
			UnassignedSequencesSpan,
			UnassignedSequencesScore.ActualSingleCopyPercentage,
			UnassignedSequencesScore.ActualDuplicatePercentage,
		)
		fmt.Println("---------------------------------------------------")
	}
}

func KmerMode(FastaFilePath string, TabexFilePath string, K int, FaiPath string, BUSCOFullTable string, NumberOfBins int, ThresholdsToTest []float32, MinKCount int, MaxKCount int) {

	// === STEP 1: Read BUSCO Full Table ===
	BUSCOstart := time.Now()
	BUSCOSequenceMap, SequenceHashMap, _, BUSCOMap := BUSCO.ReadBUSCOFullTable(BUSCOFullTable)
	TotalBUSCONum := len(BUSCOMap)
	elapsed := time.Since(BUSCOstart)
	fmt.Printf("Time elapsed to process BUSCO files: %s\n", elapsed)

	// === STEP 2: Read sequence lengths ===
	faistart := time.Now()
	SequenceLengthMap := Fai.ReadFaiFile(FaiPath, SequenceHashMap)
	faielapsed := time.Since(faistart)
	fmt.Printf("Time elapsed to process Samtools fai file: %s\n", faielapsed)

	// === STEP 3: Build Kmer Map ===
	Tabexstart := time.Now()
	TabexMap := FASTK.ReadTabexFile(TabexFilePath, K, MinKCount, MaxKCount)
	Tabexelapsed := time.Since(Tabexstart)
	fmt.Printf("Time elapsed to process Tabex file: %s\n", Tabexelapsed)

	LocateKmersstart := time.Now()
	FilteredKmerList := Kmers.KmeriseWithWorkerPool(FastaFilePath, K, TabexMap, SequenceHashMap, 8)
	Locatekmerselapsed := time.Since(LocateKmersstart)
	fmt.Printf("Time elapsed to locate kmers in genome: %s\n", Locatekmerselapsed)

	KmerMapsGenerationStart := time.Now()
	KmerMap := Kmers.ConvertKmerListToKmerMap(FilteredKmerList)
	TotalKmerCount := len(KmerMap)
	KmerSequenceMap := Kmers.GenerateSequenceMap(KmerMap)
	KmerMapsGenerationelapsed := time.Since(KmerMapsGenerationStart)
	fmt.Printf("Time elapsed to generate KmerMap and KmerSequenceMap: %s\n", KmerMapsGenerationelapsed)

	// === STEP 4: Construct graph from k-mers ===
	graphtimestart := time.Now()
	BufferGraph, graph := Graph.GenerateGraphInParallel(KmerMap, 8)
	graphelapsed := time.Since(graphtimestart)
	fmt.Printf("Time elapsed to build graph from kmers: %s\n", graphelapsed)

	graphpartitiontimestart := time.Now()
	SortedSequences := Graph.SortAllSequencesbyConnectedness(graph)

	fmt.Println("Threshold\tBin\tSeqs\tSpan\tS(%)\tD(%)")

	for _, Threshold := range ThresholdsToTest {
		// === STEP 5: Bin sequences using graph ===
		Bins, UnassignedSequences := Graph.PartitionGraph(
			BufferGraph,
			graph,
			KmerSequenceMap,
			SortedSequences,
			NumberOfBins,
			Threshold,
			TotalKmerCount,
			SequenceLengthMap,
		)

		// === STEP 6: Score and report each bin using BUSCOs ===
		for i, Bin := range Bins {
			var BinBUSCOs []uint64
			for _, Sequence := range Bin.Sequences {
				BinBUSCOs = append(BinBUSCOs, BUSCOSequenceMap[Sequence].Markers...)
			}

			BinScore := Sequence.ScoreMarkers(BinBUSCOs, TotalBUSCONum)

			fmt.Printf("%.1f\t%d\t%d\t%d\t%.2f\t%.2f\n",
				Threshold,
				i+1,
				Bin.NumofSequences,
				Bin.AssemblySpan,
				BinScore.ActualSingleCopyPercentage,
				BinScore.ActualDuplicatePercentage,
			)
		}

		// === STEP 7: Score unassigned sequences ===
		var UnassignedSequencesBUSCOs []uint64
		for _, Sequence := range UnassignedSequences {
			UnassignedSequencesBUSCOs = append(UnassignedSequencesBUSCOs, BUSCOSequenceMap[Sequence].Markers...)
		}
		UnassignedSequencesScore := Sequence.ScoreMarkers(UnassignedSequencesBUSCOs, TotalBUSCONum)
		UnassignedSequencesSpan := Fai.MeasureAssemblySpan(UnassignedSequences, SequenceLengthMap)

		fmt.Printf("Unassigned\t%d\t%d\t%.2f\t%.2f\n",
			len(UnassignedSequences),
			UnassignedSequencesSpan,
			UnassignedSequencesScore.ActualSingleCopyPercentage,
			UnassignedSequencesScore.ActualDuplicatePercentage,
		)
		fmt.Println("---------------------------------------------------")
	}
	graphpartitionelapsed := time.Since(graphpartitiontimestart)
	fmt.Printf("Time elapsed to partition graph and run sub-processes: %s\t", graphpartitionelapsed)
}

func main() {

	FastaPath := flag.String("Fasta", "", "Unphased fasta file")
	TabexPath := flag.String("Tabex", "", "Tabex -A -t output for the unphased fasta")
	FaiPath := flag.String("Fai", "", "Samtools index .fai file for the unphased fasta")
	K := flag.Int("K", 0, "Kmer length used when running FastK")
	BinCount := flag.Int("Ploidy", 0, "Ploidy, or number of bins you want to partition your genome into")
	Min := flag.Int("MinKCount", 2, "Minimum genome-wide count of kmers used in phasing, default: 2")
	Max := flag.Int("MaxKCount", -1, "Maximum genome-wide count of kmers used in phasing, default: Ploidy")
	BuscoTable := flag.String("BUSCOFullTable", "", "BUSCO full_table.tsv file for the unphased fasta")
	flag.Parse()

	// Check required flags
	if *FastaPath == "" || *TabexPath == "" || *FaiPath == "" || *K == 0 || *BinCount == 0 {
		fmt.Fprintln(os.Stderr, "Error: Required flags -Fasta, -Tabex, -Fai, -K, and -Ploidy must be provided as non-zero/non-empty")
		flag.Usage()
		os.Exit(1)
	}

	if *Max < *Min {
		fmt.Println("MaxKCount must be bigger than or equal to MinKCount")
		os.Exit(1)
	}

	// If MaxKCount wasn't provided, default to Ploidy
	if *Max == -1 {
		*Max = *BinCount
	}

	Thresholds := []float32{5, 10, 15}

	KmerMode(*FastaPath, *TabexPath, *K, *FaiPath, *BuscoTable, *BinCount, Thresholds, *Min, *Max)

}

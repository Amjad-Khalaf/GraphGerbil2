package main

import (
	"GraphGerbil2/BUSCO"
	"GraphGerbil2/Fai"
	"GraphGerbil2/Graph"
	"GraphGerbil2/Sequence"
	"GraphGerbil2/Syncmers"
	"flag"
	"fmt"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
)


func BUSCOMode(FaiPath string, BUSCOFullTable string, NumberOfBins int, Thresholds []float32, NumOfWorkers int) {

	start := time.Now()
	SequenceLengthMap, SequenceHashMap := Fai.ReadFaiFile(FaiPath)
	elapsed := time.Since(start)
	fmt.Printf("Time elapsed to process Samtools fai file: %s\n", elapsed)
	
	start = time.Now()
	BUSCOSequenceMap, _, BUSCOMap, TotalBUSCONum := BUSCO.ReadBUSCOFullTable(BUSCOFullTable, SequenceHashMap)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to process BUSCO files: %s\n", elapsed)

	fmt.Printf("%d sequences in total, %d of those have BUSCOs on them \n", len(SequenceLengthMap), len(BUSCOSequenceMap))

	start = time.Now()
	BufferGraph, graph := Graph.GenerateGraphInParallel(BUSCOMap, NumOfWorkers)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to build graph: %s\n", elapsed)

	graphpartitiontimestart := time.Now()
	SortedSequences := Graph.SortAllSequencesbyConnectedness(graph)

	fmt.Println("---------------------------------------------------")
	fmt.Printf("%-10s %-5s %-5s %-10s %-6s %-6s %-6s\n", "Threshold", "Bin", "Seqs", "Span", "C(%)", "S(%)", "D(%)")
	fmt.Println("---------------------------------------------------")


	for _, Threshold := range Thresholds {
		Bins, UnassignedSequences := Graph.PartitionGraph(
			BufferGraph,
			graph,
			BUSCOSequenceMap,
			SortedSequences,
			NumberOfBins,
			Threshold,
			TotalBUSCONum,
			SequenceLengthMap,
		)

		for i, Bin := range Bins {
			fmt.Printf("%-10.1f %-5d %-5d %-10d %-6.2f %-6.2f %-6.2f\n",
				Threshold,
				i+1,
				Bin.NumofSequences,
				Bin.AssemblySpan,
				(Bin.Score.ActualSingleCopyPercentage + Bin.Score.ActualDuplicatePercentage),
				Bin.Score.ActualSingleCopyPercentage,
				Bin.Score.ActualDuplicatePercentage,
			)
		}

		var UnassignedSequencesBUSCOs []uint64
		for _, Sequence := range UnassignedSequences {
			UnassignedSequencesBUSCOs = append(UnassignedSequencesBUSCOs, BUSCOSequenceMap[Sequence].Markers...)
		}
		UnassignedSequencesScore := Sequence.ScoreMarkers(UnassignedSequencesBUSCOs, TotalBUSCONum)
		UnassignedSequencesSpan := Fai.MeasureAssemblySpan(UnassignedSequences, SequenceLengthMap)

		fmt.Printf("%-10s %-5s %-5d %-10d %-6.2f %-6.2f %-6.2f\n",
			"Unassigned",
			"",
			len(UnassignedSequences),
			UnassignedSequencesSpan,
			(UnassignedSequencesScore.ActualSingleCopyPercentage + UnassignedSequencesScore.ActualDuplicatePercentage),
			UnassignedSequencesScore.ActualSingleCopyPercentage,
			UnassignedSequencesScore.ActualDuplicatePercentage,
		)
		fmt.Println("---------------------------------------------------")
	}
	graphpartitionelapsed := time.Since(graphpartitiontimestart)
	fmt.Printf("Time elapsed to partition graph and run sub-processes: %s\t", graphpartitionelapsed)
}

func SyncmerMode(FastaFilePath string, K int, S int, T int, FaiPath string, BUSCOFullTable string, NumberOfBins int, Thresholds []float32, MinKCount int, MaxKCount int, cores int, NumToSubsample int) {

	start := time.Now()
	SequenceLengthMap, SequenceHashMap := Fai.ReadFaiFile(FaiPath)
	elapsed := time.Since(start)
	fmt.Printf("Time elapsed to process Samtools fai file: %s\n", elapsed)

	start = time.Now()
	BUSCOSequenceMap, _, _, TotalBUSCONum := BUSCO.ReadBUSCOFullTable(BUSCOFullTable, SequenceHashMap)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to process BUSCO files: %s\n", elapsed)

	fmt.Printf("%d sequences in total, %d of those have BUSCOs on them \n", len(SequenceLengthMap), len(BUSCOSequenceMap))

	start = time.Now()
	SyncmerList := Syncmers.SyncmeriseWithWorkerPool(FastaFilePath, K, S, T, SequenceHashMap, cores)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to generate Syncmers for the whole genome: %s\n", elapsed)

	start = time.Now()
	SyncmerMap := Syncmers.ConvertSyncmerListToKSyncmerMapIgnoringInternalDuplication(SyncmerList)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to generate SyncmerMap from SyncmerList: %s\n", elapsed)

	start = time.Now()
	FilteredSyncmerMap := Syncmers.FilterSyncmerMap(SyncmerMap, MinKCount, MaxKCount)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to filter syncmers: %s\n", elapsed)

	TotalSyncmerCount := len(FilteredSyncmerMap)

	start = time.Now()
	SyncmerSequenceMap := Syncmers.GenerateSequenceMap(FilteredSyncmerMap)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to generate SyncmerSequenceMap from Filtered SyncmerMap: %s\n", elapsed)
	fmt.Printf("%d sequences in total are represented by %d valid syncmers \n", len(SyncmerSequenceMap), TotalSyncmerCount)

	if NumToSubsample != -1 {
		start = time.Now()
		SampledSyncmerMap := Syncmers.SampleSyncmers(FilteredSyncmerMap, NumToSubsample)
		elapsed = time.Since(start)
		fmt.Printf("Time elapsed to sample %d syncmers: %s\n", NumToSubsample, elapsed)
		SyncmerSequenceMap = Syncmers.GenerateSequenceMap(SampledSyncmerMap)
		fmt.Printf("%d sequences in total are represented by sampled syncmers \n", len(SyncmerSequenceMap))
		FilteredSyncmerMap = SampledSyncmerMap
	}

	start = time.Now()
	BufferGraph, graph := Graph.GenerateGraphInParallel(FilteredSyncmerMap, cores)
	elapsed = time.Since(start)
	fmt.Printf("Time elapsed to build graph: %s\n", elapsed)

	graphpartitiontimestart := time.Now()
	SortedSequences := Graph.SortAllSequencesbyConnectedness(graph)

	fmt.Println("---------------------------------------------------")
	fmt.Printf("%-10s %-5s %-5s %-10s %-6s %-6s %-6s\n", "Threshold", "Bin", "Seqs", "Span", "C(%)", "S(%)", "D(%)")
	fmt.Println("---------------------------------------------------")

	for _, Threshold := range Thresholds {
		Bins, UnassignedSequences := Graph.PartitionGraphInParallel(
			BufferGraph,
			graph,
			SyncmerSequenceMap,
			SortedSequences,
			NumberOfBins,
			Threshold,
			TotalSyncmerCount,
			SequenceLengthMap,
		)

		for i, Bin := range Bins {
			var BinBUSCOs []uint64
			for _, Sequence := range Bin.Sequences {
				BinBUSCOs = append(BinBUSCOs, BUSCOSequenceMap[Sequence].Markers...)
			}

			BinScore := Sequence.ScoreMarkers(BinBUSCOs, TotalBUSCONum)

			fmt.Printf("%-10.1f %-5d %-5d %-10d %-6.2f %-6.2f %-6.2f\n",
				Threshold,
				i+1,
				Bin.NumofSequences,
				Bin.AssemblySpan,
				(BinScore.ActualSingleCopyPercentage + BinScore.ActualDuplicatePercentage),
				BinScore.ActualSingleCopyPercentage,
				BinScore.ActualDuplicatePercentage,
			)
		}

		var UnassignedSequencesBUSCOs []uint64
		for _, Sequence := range UnassignedSequences {
			UnassignedSequencesBUSCOs = append(UnassignedSequencesBUSCOs, BUSCOSequenceMap[Sequence].Markers...)
		}
		UnassignedSequencesScore := Sequence.ScoreMarkers(UnassignedSequencesBUSCOs, TotalBUSCONum)
		UnassignedSequencesSpan := Fai.MeasureAssemblySpan(UnassignedSequences, SequenceLengthMap)

		fmt.Printf("%-10s %-5s %-5d %-10d %-6.2f %-6.2f %-6.2f\n",
			"Unassigned",
			"",
			len(UnassignedSequences),
			UnassignedSequencesSpan,
			(UnassignedSequencesScore.ActualSingleCopyPercentage + UnassignedSequencesScore.ActualDuplicatePercentage),
			UnassignedSequencesScore.ActualSingleCopyPercentage,
			UnassignedSequencesScore.ActualDuplicatePercentage,
		)
		fmt.Println("---------------------------------------------------")
	}
	graphpartitionelapsed := time.Since(graphpartitiontimestart)
	fmt.Printf("Time elapsed to partition graph and run sub-processes: %s\t", graphpartitionelapsed)
}

func main() {
	mode := flag.String("Mode", "", "Phasing mode to use: 'syncmers' or 'busco'")
	FastaPath := flag.String("Fasta", "", "(non-BUSCO modes only) Unphased fasta file")
	FaiPath := flag.String("Fai", "", "Samtools index .fai file for the unphased fasta")
	K := flag.Int("K", 31, "(non-BUSCO modes only) Kmer length, default 31")
	S := flag.Int("S", 11, "(non-BUSCO modes only) Syncmer S, default 11")
	T := flag.Int("T", 1, "(non-BUSCO modes only) Syncmer T, default 1")
	Ploidy := flag.Int("Ploidy", 0, "Ploidy, or number of bins you want to partition your genome into")
	Min := flag.Int("Min", 2, "(non-BUSCO modes only) Minimum genome-wide count of markers used in phasing, default: 2")
	Max := flag.Int("Max", -1, "(non-BUSCO modes only) Maximum genome-wide count of markers used in phasing, default: Ploidy")
	BuscoTable := flag.String("BuscoFullTable", "", "BUSCO full_table.tsv file for the unphased fasta")
	ThresholdString := flag.String("DuplicationThresholds", "1,10,20,30,40,50", "Starting thresholds to test, please provide as comma-separated list, e.g.: 1,2,3")
	NumOfWorkers := flag.Int("cores", -1, "Number of cores you want the programme to use, default is all available")
	NumToSubsample := flag.Int("R", -1, "(non-BUSCO modes only) Number of markers to subsample when building the graph, default is to use all markers available")
	flag.Parse()

	if *mode != "syncmers" && *mode != "busco" {
		fmt.Fprintf(os.Stderr, "ERROR: Please pick a valid mode\n")
		flag.Usage()
		os.Exit(1)
	}

	missing := func(name string) {
		fmt.Fprintf(os.Stderr, "ERROR: -%s is required\n", name)
		flag.Usage()
		os.Exit(1)
	}
	if *mode == "" {
		missing("Mode")
	}
	if *FastaPath == ""  && *mode == "syncmers"{
		missing("Fasta")
	}
	if *FaiPath == "" {
		missing("Fai")
	}
	if *Ploidy == 0 {
		missing("Ploidy")
	}
	if *BuscoTable == "" {
		missing("BuscoFullTable")
	}

	cores := *NumOfWorkers
	if cores <= 0 {
		cores = runtime.NumCPU()
	}
	maxKCount := *Max
	if maxKCount < 0 {
		maxKCount = *Ploidy
	}

	var thresholds []float32
	for _, tok := range strings.Split(*ThresholdString, ",") {
		if strings.TrimSpace(tok) == "" {
			continue
		}
		val, err := strconv.ParseFloat(strings.TrimSpace(tok), 32)
		if err != nil {
			fmt.Fprintf(os.Stderr, "ERROR: cannot parse threshold %q: %v\n", tok, err)
			os.Exit(1)
		}
		thresholds = append(thresholds, float32(val))
	}
	if len(thresholds) == 0 {
		fmt.Fprintln(os.Stderr, "ERROR: no valid thresholds supplied with -DuplicationThresholds")
		os.Exit(1)
	}


	if *mode == "syncmers" {
		fmt.Printf("\nLaunching SyncmerMode on %s with %d cores …\n", *FastaPath, cores)
		SyncmerMode(
			*FastaPath,  // FastaFilePath
			*K,          // K
			*S,          // S
			*T,          // T
			*FaiPath,    // FaiPath
			*BuscoTable, // BUSCOFullTable
			*Ploidy,     // NumberOfBins
			thresholds,  // Thresholds
			*Min,        // MinKCount
			maxKCount,   // MaxKCount
			cores,       // cores
			*NumToSubsample,
		)
	}

	if *mode == "busco" {
		fmt.Printf("\nLaunching BUSCOmode on %s with %d cores …\n", *FastaPath, cores)
		BUSCOMode(
			*FaiPath,    // FaiPath
			*BuscoTable, // BUSCOFullTable
			*Ploidy,     // NumberOfBins
			thresholds,  // Thresholds
			cores,       // cores
		)
	}
}

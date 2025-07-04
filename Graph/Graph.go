package Graph

import (
	"GraphGerbil2/Fai"
	"GraphGerbil2/Sequence"
	"slices"
	"sort"
	"math"
	"sync"
)

type Edge struct {
	To     uint64 // Where edge leads to
	Weight int    // Number of shared Markers
}

type Candidate struct {
	Name                uint64
	TotalMarkerNumAdded int
}

type ScoredSequence struct {
	Name          uint64
	Score         Sequence.Score
	JoinedMarkers []uint64
}

type GenomeAssemblyBin struct {
	Sequences      []uint64
	NumofSequences int
	AssemblySpan   uint64
	Markers        []uint64
	Score          Sequence.Score
}

func GenerateGraph(MarkerMap map[uint64][]uint64) (map[uint64]map[uint64]int, map[uint64][]Edge) {
	BufferGraph := make(map[uint64]map[uint64]int) // By storing edges as a map, with Sequences as keys, look up should be faster than checking a list of edges by iteration
	for _, SequenceList := range MarkerMap {
		for i := 0; i < len(SequenceList); i++ {
			for j := i + 1; j < len(SequenceList); j++ {
				SequenceA := SequenceList[i]
				SequenceB := SequenceList[j]

				// Initialise sub-maps
				if BufferGraph[SequenceA] == nil {
					BufferGraph[SequenceA] = make(map[uint64]int)
				}
				if BufferGraph[SequenceB] == nil {
					BufferGraph[SequenceB] = make(map[uint64]int)
				}

				// Increment for shared Markers
				BufferGraph[SequenceA][SequenceB]++
				BufferGraph[SequenceB][SequenceA]++
			}
		}
	}
	Graph := make(map[uint64][]Edge)
	for SequenceA, AdjacencyMaps := range BufferGraph {
		for SequenceB, SharedMarkerNum := range AdjacencyMaps {
			Info := Edge{To: SequenceB, Weight: SharedMarkerNum}
			Graph[SequenceA] = append(Graph[SequenceA], Info)
		}
	}
	return BufferGraph, Graph
}

func Worker(SequenceList []uint64) map[uint64]map[uint64]int {
	LocalGraph := make(map[uint64]map[uint64]int)
	for i := 0; i < len(SequenceList); i++ {
		for j := i + 1; j < len(SequenceList); j++ {
			SequenceA := SequenceList[i]
			SequenceB := SequenceList[j]

			// Initialise sub-maps
			if LocalGraph[SequenceA] == nil {
				LocalGraph[SequenceA] = make(map[uint64]int)
			}
			if LocalGraph[SequenceB] == nil {
				LocalGraph[SequenceB] = make(map[uint64]int)
			}

			// Increment for shared Markers
			LocalGraph[SequenceA][SequenceB]++
			LocalGraph[SequenceB][SequenceA]++
		}
	}
	return LocalGraph
}

func GenerateGraphInParallel(MarkerMap map[uint64][]uint64, NumOfWorkers int) (map[uint64]map[uint64]int, map[uint64][]Edge) {

	// One unit of work is to take a []uint64 and generate all pairs and track co-ocurrences

	// Each "Job" is a unit of work received by a worker
	type Job struct {
		Marker       uint64
		SequenceList []uint64
	}

	// Each result is what the worker returns with (i.e. a local graph)
	// The main thread should eventually merge all these into the global BufferGraph
	type Result struct {
		LocalGraph map[uint64]map[uint64]int
	}

	// Define channels
	Jobs := make(chan Job, len(MarkerMap))
	Results := make(chan Result, NumOfWorkers)

	// Launch workers and send them jobs
	for i := 0; i < NumOfWorkers; i++ {
		go func() {
			for Job := range Jobs {
				localGraph := Worker(Job.SequenceList)
				Results <- Result{LocalGraph: localGraph}
			}
		}()
	}

	// Send all jobs
	for Marker, SequenceList := range MarkerMap {
		Jobs <- Job{Marker: Marker, SequenceList: SequenceList}
	}

	// Close Jobs channel
	close(Jobs)

	// Collect all results and merge them into global BufferGraph
	BufferGraph := make(map[uint64]map[uint64]int)
	for i := 0; i < len(MarkerMap); i++ {
		MergedResult := <-Results

		for SequenceA, AdjacencyMap := range MergedResult.LocalGraph {
			if BufferGraph[SequenceA] == nil {
				BufferGraph[SequenceA] = make(map[uint64]int)
			}
			for SequenceB, Count := range AdjacencyMap {
				BufferGraph[SequenceA][SequenceB] += Count
			}
		}
	}

	// Convert global buffer graph to graph
	Graph := make(map[uint64][]Edge)
	for SequenceA, AdjacencyMaps := range BufferGraph {
		for SequenceB, SharedMarkerNum := range AdjacencyMaps {
			Info := Edge{To: SequenceB, Weight: SharedMarkerNum}
			Graph[SequenceA] = append(Graph[SequenceA], Info)
		}
	}
	return BufferGraph, Graph
}

func SortAllSequencesbyConnectedness(Graph map[uint64][]Edge) []uint64 {
	// Sorts sequence list by connectivity, in descending order
	var SequenceList []uint64
	var NumOfConnections []int
	for Sequence, EdgeList := range Graph {
		SequenceList = append(SequenceList, Sequence)
		NumOfConnections = append(NumOfConnections, len(EdgeList))
	}

	var SortedSequenceList []uint64
	SortedNumOfConnections := slices.Clone(NumOfConnections)
	slices.Sort(SortedNumOfConnections)
	slices.Reverse(SortedNumOfConnections)
	var OriginalIndex int
	var UsedIndeces []int
	for _, Length := range SortedNumOfConnections {
		for index, value := range NumOfConnections {
			if Length == value && !slices.Contains(UsedIndeces, index) {
				OriginalIndex = index
				UsedIndeces = append(UsedIndeces, index)
				break
			}
		}
		SortedSequenceList = append(SortedSequenceList, SequenceList[OriginalIndex])
	}
	return SortedSequenceList
}

func SortSequencesSubsetbyConnectedness(Graph map[uint64][]Edge, Subset []uint64) []uint64 {
	Scores := make(map[uint64]int)
	for _, Sequence := range Subset {
		Scores[Sequence] = len(Graph[Sequence])
	}
	sort.Slice(Subset, func(i, j int) bool {
		return Scores[Subset[i]] > Scores[Subset[j]] // descending
	})
	return Subset // Sorted
}

func PartitionGraph(BufferGraph map[uint64]map[uint64]int, Graph map[uint64][]Edge, SequenceMap map[uint64]Sequence.Sequence, SequenceList []uint64, BinNum int, RelativeDuplicationThreshold float32, TotalMarkerNum int, SequenceLengthMap map[uint64]uint64) ([]GenomeAssemblyBin, []uint64) {
	Assigned := make(map[uint64]bool)
	var GenomeAssemblyBins []GenomeAssemblyBin

	for BinIndex := 0; BinIndex < BinNum; BinIndex++ {
		// Get unassigned sequences
		var UnassignedSequences []uint64
		for _, seq := range SequenceList {
			if !Assigned[seq] {
				UnassignedSequences = append(UnassignedSequences, seq)
			}
		}
		if len(UnassignedSequences) == 0 {
			break
		}

		// Initialize bin with the most connected sequence
		Sorted := SortSequencesSubsetbyConnectedness(Graph, UnassignedSequences)
		seed := Sorted[0]
		binSequences := []uint64{seed}
		binMarkers := slices.Clone(SequenceMap[seed].Markers) // Avoid mutation
		Assigned[seed] = true

		// Track Marker counts
		MarkerCounts := make(map[uint64]int)
		for _, Marker := range binMarkers {
			MarkerCounts[Marker]++
		}

		// Counters for the bin
		singleCopy := 0
		duplicated := 0
		for _, count := range MarkerCounts {
			if count == 1 {
				singleCopy++
			} else if count > 1 {
				duplicated++
			}
		}

		// Candidate processing
		candidates := Sorted[1:]

		for len(candidates) > 0 {
			var bestCandidate uint64 = 100000000
			bestNewDuplicates := int(^uint(0) >> 1) // Max int

			for _, candidate := range candidates {
				if Assigned[candidate] {
					continue
				}

				// Simulate adding candidate
				newDuplicates := 0
				for _, Marker := range SequenceMap[candidate].Markers {
					currentCount := MarkerCounts[Marker]
					switch currentCount {
					case 1:
						newDuplicates++ // Becomes duplicated
					case 0:
						// New Marker, no duplicate
					default:
						// Already duplicated, no change
					}
				}

				// Skip if exceeds threshold
				newDuplicated := duplicated
				if newDuplicates > 0 {
					newDuplicated += newDuplicates
				}
				newTotal := singleCopy + duplicated + len(SequenceMap[candidate].Markers) - newDuplicates
				newDupPercent := (float32(newDuplicated) / float32(newTotal) * 100)

				if newDupPercent <= RelativeDuplicationThreshold && newDuplicates < bestNewDuplicates {
					bestNewDuplicates = newDuplicates
					bestCandidate = candidate
				}
			}

			if bestCandidate == 100000000 {
				break // No valid candidate
			}

			// Add best candidate to bin
			binSequences = append(binSequences, bestCandidate)
			Assigned[bestCandidate] = true
			for _, Marker := range SequenceMap[bestCandidate].Markers {
				oldCount := MarkerCounts[Marker]
				newCount := oldCount + 1
				MarkerCounts[Marker] = newCount

				if oldCount == 1 {
					duplicated++ // Now duplicated
					singleCopy--
				} else if oldCount == 0 {
					singleCopy++ // New single-copy
				}
				// (No change if already duplicated)
			}

			// Remove candidate from list
			newCandidates := make([]uint64, 0, len(candidates)-1)
			for _, c := range candidates {
				if c != bestCandidate {
					newCandidates = append(newCandidates, c)
				}
			}
			candidates = newCandidates
		}

		// Build final Marker list and score
		allMarkers := make([]uint64, 0, len(MarkerCounts)*2) // Estimate size
		for Marker, count := range MarkerCounts {
			for i := 0; i < count; i++ {
				allMarkers = append(allMarkers, Marker)
			}
		}
		finalScore := Sequence.ScoreMarkers(allMarkers, TotalMarkerNum)

		bin := GenomeAssemblyBin{
			Sequences:      binSequences,
			NumofSequences: len(binSequences),
			AssemblySpan:   Fai.MeasureAssemblySpan(binSequences, SequenceLengthMap),
			Markers:        allMarkers,
			Score:          finalScore,
		}
		GenomeAssemblyBins = append(GenomeAssemblyBins, bin)
	}

	var FinalUnassignedSequences []uint64
	for _, Sequence := range SequenceList {
		if !Assigned[Sequence] {
			FinalUnassignedSequences = append(FinalUnassignedSequences, Sequence)
		}
	}
	return GenomeAssemblyBins, FinalUnassignedSequences
}

func PartitionGraphInParallel(
	BufferGraph map[uint64]map[uint64]int,         // (unused here; kept for API stability)
	Graph map[uint64][]Edge,
	SequenceMap map[uint64]Sequence.Sequence,
	SequenceList []uint64,
	BinNum int,
	RelativeDuplicationThreshold float32,
	TotalMarkerNum int,
	SequenceLengthMap map[uint64]uint64,
) ([]GenomeAssemblyBin, []uint64) {

	const minParallelCandidates = 4 // don’t incur goroutine overhead for tiny sets

	Assigned := make(map[uint64]bool)
	var GenomeAssemblyBins []GenomeAssemblyBin

	for BinIndex := 0; BinIndex < BinNum; BinIndex++ {

		// 1. gather unassigned sequences -------------------------------------
		var Unassigned []uint64
		for _, seq := range SequenceList {
			if !Assigned[seq] {
				Unassigned = append(Unassigned, seq)
			}
		}
		if len(Unassigned) == 0 {
			break
		}

		// 2. seed the bin with the most‑connected sequence --------------------
		sorted := SortSequencesSubsetbyConnectedness(Graph, Unassigned)
		seed := sorted[0]

		binSeqs   := []uint64{seed}
		binMarks  := slices.Clone(SequenceMap[seed].Markers) // defensive copy
		Assigned[seed] = true

		// 3. marker accounting -------------------------------------------------
		markCounts := make(map[uint64]int, len(binMarks))
		for _, m := range binMarks {
			markCounts[m]++
		}

		singleCopy, duplicated := 0, 0
		for _, c := range markCounts {
			if c == 1 {
				singleCopy++
			} else if c > 1 {
				duplicated++
			}
		}

		// 4. fill the bin ------------------------------------------------------
		candidates := sorted[1:]

		for len(candidates) > 0 {

			bestCand      := uint64(0)       // 0 => none found
			bestNewDup    := math.MaxInt
			bestCandIndex := -1

			// -------- parallel candidate evaluation --------------------------
			if len(candidates) >= minParallelCandidates {
				type candEval struct {
					id            uint64
					newDup        int
					newDupPercent float32
					index         int
				}

				results := make(chan candEval, len(candidates))
				var wg sync.WaitGroup

				// immutable snapshots for goroutines
				snapCounts   := markCounts
				curDup       := duplicated
				curSingle    := singleCopy

				for idx, cand := range candidates {
					if Assigned[cand] { // should not happen, but cheap guard
						continue
					}
					wg.Add(1)
					go func(i int, c uint64) {
						defer wg.Done()

						newDup := 0
						for _, m := range SequenceMap[c].Markers {
							if snapCounts[m] == 1 {
								newDup++ // will flip single→dup
							}
						}

						newDuplicated := curDup + newDup
						newTotal := curSingle + curDup +
							len(SequenceMap[c].Markers) - newDup
						dupPct := float32(newDuplicated) / float32(newTotal) * 100

						results <- candEval{c, newDup, dupPct, i}
					}(idx, cand)
				}

				wg.Wait()
				close(results)

				for r := range results {
					if r.newDupPercent > RelativeDuplicationThreshold {
						continue
					}
					if r.newDup < bestNewDup ||
						(r.newDup == bestNewDup && r.index < bestCandIndex) {
						bestCand, bestNewDup, bestCandIndex =
							r.id, r.newDup, r.index
					}
				}
			}

			// -------- fallback: tiny candidate set, serial scan --------------
			if bestCand == 0 { // either we didn’t parallelise or none passed
				for idx, cand := range candidates {
					if Assigned[cand] {
						continue
					}

					newDup := 0
					for _, m := range SequenceMap[cand].Markers {
						if markCounts[m] == 1 {
							newDup++
						}
					}

					newDuplicated := duplicated + newDup
					newTotal := singleCopy + duplicated +
						len(SequenceMap[cand].Markers) - newDup
					dupPct := float32(newDuplicated) / float32(newTotal) * 100
					if dupPct > RelativeDuplicationThreshold {
						continue
					}

					if newDup < bestNewDup ||
						(newDup == bestNewDup && (bestCandIndex == -1 || idx < bestCandIndex)) {
						bestCand, bestNewDup, bestCandIndex =
							cand, newDup, idx
					}
				}
			}

			// no candidate satisfies the threshold → finish this bin ----------
			if bestCand == 0 {
				break
			}

			// -------- mutate state with the chosen candidate -----------------
			binSeqs = append(binSeqs, bestCand)
			Assigned[bestCand] = true

			for _, m := range SequenceMap[bestCand].Markers {
				switch markCounts[m] {
				case 0:
					singleCopy++
				case 1:
					singleCopy--
					duplicated++
				}
				markCounts[m]++
			}

			// remove from candidate list without extra alloc
			candidates = append(candidates[:bestCandIndex],
				candidates[bestCandIndex+1:]...)
		}

		// 5. build bin object --------------------------------------------------
		var allMarkers []uint64
		for m, c := range markCounts {
			for i := 0; i < c; i++ {
				allMarkers = append(allMarkers, m)
			}
		}

		score := Sequence.ScoreMarkers(allMarkers, TotalMarkerNum)

		GenomeAssemblyBins = append(GenomeAssemblyBins, GenomeAssemblyBin{
			Sequences:      binSeqs,
			NumofSequences: len(binSeqs),
			AssemblySpan:   Fai.MeasureAssemblySpan(binSeqs, SequenceLengthMap),
			Markers:        allMarkers,
			Score:          score,
		})
	}

	// 6. collect any sequences that never fit into a bin ----------------------
	var unassigned []uint64
	for _, seq := range SequenceList {
		if !Assigned[seq] {
			unassigned = append(unassigned, seq)
		}
	}

	return GenomeAssemblyBins, unassigned
}

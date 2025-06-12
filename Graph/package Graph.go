package Graph

import (
	"GraphGerbil/Sequence"
	"GraphGerbil/fai"
	"slices"
	"sort"
)

type Edge struct {
	To     string // Where edge leads to
	Weight int    // Number of shared Markers
}

type Candidate struct {
	Name                string
	TotalMarkerNumAdded int
}

type ScoredSequence struct {
	Name          string
	Score         Sequence.Score
	JoinedMarkers []string
}

type GenomeAssemblyBin struct {
	Sequences      []string
	NumofSequences int
	AssemblySpan   int
	Markers        []string
	Score          Sequence.Score
}

func GenerateGraph(MarkerMap map[string][]string) (map[string]map[string]int, map[string][]Edge) {
	// MarkerMap stores Sequence names as keys, and a slice of their Markers as values
	BufferGraph := make(map[string]map[string]int) // By storing edges as a map, with Sequences as keys, look up should be faster than checking a list of edges by iteration
	for _, SequenceList := range MarkerMap {
		for i := 0; i < len(SequenceList); i++ {
			for j := i + 1; j < len(SequenceList); j++ {
				SequenceA := SequenceList[i]
				SequenceB := SequenceList[j]

				// Initialise sub-maps
				if BufferGraph[SequenceA] == nil {
					BufferGraph[SequenceA] = make(map[string]int)
				}
				if BufferGraph[SequenceB] == nil {
					BufferGraph[SequenceB] = make(map[string]int)
				}

				// Increment for shared Markers
				BufferGraph[SequenceA][SequenceB]++
				BufferGraph[SequenceB][SequenceA]++
			}
		}
	}
	Graph := make(map[string][]Edge)
	for SequenceA, AdjacencyMaps := range BufferGraph {
		for SequenceB, SharedMarkerNum := range AdjacencyMaps {
			Info := Edge{To: SequenceB, Weight: SharedMarkerNum}
			Graph[SequenceA] = append(Graph[SequenceA], Info)
		}
	}
	return BufferGraph, Graph
}

func SortAllSequencesbyConnectedness(Graph map[string][]Edge) []string {
	// Sorts sequence list by connectivity, in descending order
	var SequenceList []string
	var NumOfConnections []int
	for Sequence, EdgeList := range Graph {
		SequenceList = append(SequenceList, Sequence)
		NumOfConnections = append(NumOfConnections, len(EdgeList))
	}

	var SortedSequenceList []string
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

func SortSequencesSubsetbyConnectedness(Graph map[string][]Edge, Subset []string) []string {
	Scores := make(map[string]int)
	for _, Sequence := range Subset {
		Scores[Sequence] = len(Graph[Sequence])
	}
	sort.Slice(Subset, func(i, j int) bool {
		return Scores[Subset[i]] > Scores[Subset[j]] // descending
	})
	return Subset // Sorted
}

func PartitionGraph(BufferGraph map[string]map[string]int, Graph map[string][]Edge, SequenceMap map[string]Sequence.Sequence, SequenceList []string, BinNum int, RelativeDuplicationThreshold float32, TotalMarkerNum int, SequenceLengthMap map[string]int) ([]GenomeAssemblyBin, []string) {
	Assigned := make(map[string]bool)
	var GenomeAssemblyBins []GenomeAssemblyBin

	for BinIndex := 0; BinIndex < BinNum; BinIndex++ {
		// Get unassigned sequences
		var UnassignedSequences []string
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
		binSequences := []string{seed}
		binMarkers := slices.Clone(SequenceMap[seed].Markers) // Avoid mutation
		Assigned[seed] = true

		// Track Marker counts for efficient updates
		MarkerCounts := make(map[string]int)
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
		//slices.Reverse(candidates)

		for len(candidates) > 0 {
			bestCandidate := ""
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

			if bestCandidate == "" {
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
			newCandidates := make([]string, 0, len(candidates)-1)
			for _, c := range candidates {
				if c != bestCandidate {
					newCandidates = append(newCandidates, c)
				}
			}
			candidates = newCandidates
		}

		// Build final Marker list and score
		allMarkers := make([]string, 0, len(MarkerCounts)*2) // Estimate size
		for Marker, count := range MarkerCounts {
			for i := 0; i < count; i++ {
				allMarkers = append(allMarkers, Marker)
			}
		}
		finalScore := Sequence.ScoreMarkers(allMarkers, TotalMarkerNum)

		bin := GenomeAssemblyBin{
			Sequences:      binSequences,
			NumofSequences: len(binSequences),
			AssemblySpan:   fai.MeasureAssemblySpan(binSequences, SequenceLengthMap),
			Markers:        allMarkers,
			Score:          finalScore,
		}
		GenomeAssemblyBins = append(GenomeAssemblyBins, bin)
	}

	var FinalUnassignedSequences []string
	for _, Sequence := range SequenceList {
		if !Assigned[Sequence] {
			FinalUnassignedSequences = append(FinalUnassignedSequences, Sequence)
		}
	}
	return GenomeAssemblyBins, FinalUnassignedSequences
}

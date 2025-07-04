package Minimisers

import (
	"GraphGerbil2/Fasta"
	"GraphGerbil2/Sequence"
)

type MinimiserType struct {
	Minimiser uint64
	Position  uint64
}

var BaseToInt = [256]uint8{
	'A': 0, 'C': 1, 'G': 2, 'T': 3,
	'a': 0, 'c': 1, 'g': 2, 't': 3,
}

const Base = 4

var Complement = map[byte]byte{
	'A': 'T',
	'C': 'G',
	'G': 'C',
	'T': 'A',
}

func MinimiseCanonically(Sequence Fasta.Fasta, K int, Window int, HeaderID uint64, MinimserList []MinimiserType) []MinimiserType {
	// Roll through sequences and test
	// Maintain a hash for forward Minimiser, and one for ReverseComplement Minimiser
	// Check if the canonical Minimiser is found in Filtered Tabex List
	// If so, store Canonical Minimiser in Minimisers slice as MinimiserType

	SequenceLength := len(Sequence.Sequence)
	var CanonicalKmersList []uint64

	var ForwardHash, ReverseHash uint64 = 0, 0
	var ValidBases int = 0
	for i := 0; i < SequenceLength; i++ {
		BaseChar := Sequence.Sequence[i]
		BaseValue := uint64(BaseToInt[BaseChar])
		if BaseValue > 3 {
			// Reset if I come across illegal characters
			ForwardHash, ReverseHash = 0, 0
			ValidBases = 0
			continue
		}

		// Update forward hash: left shift by 2 and add new base
		ForwardHash = ((ForwardHash << 2) | BaseValue) & ((1 << (2 * K)) - 1)

		// Update reverse hash: right shift and add new complement base to left
		Comp := uint64(BaseToInt[Complement[BaseChar]])
		ReverseHash = (ReverseHash >> 2) | (Comp << (2 * (K - 1)))
		ValidBases++
		if ValidBases >= K {
			Canonical := ForwardHash
			if ReverseHash < ForwardHash {
				Canonical = ReverseHash
			}
			CanonicalKmersList = append(CanonicalKmersList, Canonical)
		}

		if len(CanonicalKmersList) == Window {
			// Find Smallest Minimiser
			SmallestMinimiser := CanonicalKmersList[0]
			for _, v := range CanonicalKmersList[1:] {
				if v < SmallestMinimiser {
					SmallestMinimiser = v
				}
			}

			// Save it as our minimser
			var Minimiser MinimiserType
			Minimiser.Minimiser = SmallestMinimiser
			Minimiser.Position = HeaderID
			MinimserList = append(MinimserList, Minimiser)

			// Shift CanonicalKmersList by 1
			CanonicalKmersList = CanonicalKmersList[1:]

		}
	}
	return MinimserList
}

func MinimiseWithWorkerPool(FastaFilePath string, K int, Window int, SequenceHashMap map[string]uint64, numWorkers int) []MinimiserType {
	var MinimiserList []MinimiserType
	FastaMap := Fasta.ReadFastaFile(FastaFilePath)

	type Job struct {
		Fasta    Fasta.Fasta
		HeaderID uint64
	}

	jobs := make(chan Job)
	results := make(chan []MinimiserType)

	// Workers
	for w := 0; w < numWorkers; w++ {
		go func() {
			for job := range jobs {
				localMinimisers := MinimiseCanonically(job.Fasta, K, Window, job.HeaderID, nil)
				results <- localMinimisers
			}
		}()
	}

	// Feeding jobs
	go func() {
		for _, fasta := range FastaMap {
			seqNum := SequenceHashMap[fasta.Header]
			jobs <- Job{Fasta: fasta, HeaderID: seqNum}
			seqNum++
		}
		close(jobs)
	}()

	// Collecting results
	totalSeqs := len(FastaMap)
	for i := 0; i < totalSeqs; i++ {
		localMinimisers := <-results
		MinimiserList = append(MinimiserList, localMinimisers...)
	}
	return MinimiserList
}

func ConvertMinimiserListToKMinimiserMap(MinimiserList []MinimiserType) map[uint64][]uint64 {
	MinimiserMap := make(map[uint64][]uint64)
	for _, Minimiser := range MinimiserList {
		MinimiserMap[Minimiser.Minimiser] = append(MinimiserMap[Minimiser.Minimiser], Minimiser.Position)
	}
	return MinimiserMap
}

func FilterMinimiserMap(MinimiserMap map[uint64][]uint64, MinMCount int, MaxMCount int) map[uint64][]uint64 {
	FilteredMinimiserMap := make(map[uint64][]uint64)
	for Minimiser, Positions := range MinimiserMap {
		Count := len(Positions)
		if Count >= MinMCount && Count <= MaxMCount {
			FilteredMinimiserMap[Minimiser] = Positions
		}
	}
	return FilteredMinimiserMap
}

func GenerateSequenceMap(MinimiserMap map[uint64][]uint64) map[uint64]Sequence.Sequence {
	SequenceMap := make(map[uint64]Sequence.Sequence)
	for Minimiser, SequenceList := range MinimiserMap {
		for _, Sequence := range SequenceList {
			info := SequenceMap[Sequence]
			info.Markers = append(info.Markers, Minimiser)
			SequenceMap[Sequence] = info
		}
	}
	for Sequence, SequenceInfo := range SequenceMap {
		info := SequenceInfo
		info.NumofMarkers = len(info.Markers)
		SequenceMap[Sequence] = info
	}
	return SequenceMap
}

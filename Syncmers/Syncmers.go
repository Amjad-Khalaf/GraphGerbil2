package Syncmers

import (
	"GraphGerbil2/Fasta"
	"GraphGerbil2/Sequence"
	"math/rand"
	"time"
	"slices"
)

type SyncmerType struct {
	Syncmer  uint64
	Position uint64
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

func SyncmeriseCanonically(Sequence Fasta.Fasta, K int, S int, T int, HeaderID uint64, SyncmerList []SyncmerType) []SyncmerType {
	SequenceLength := len(Sequence.Sequence)

	var SmerForwardHash, SmerReverseHash uint64
	SmerBuffer := make([]uint64, K-S+1) // Ring buffer size is k - s + 1
	SmerPosition := 0

	var ForwardHash, ReverseHash uint64 = 0, 0
	var ValidBases int = 0

	for i := 0; i < SequenceLength; i++ {
		BaseChar := Sequence.Sequence[i]
		BaseValue := uint64(BaseToInt[BaseChar])
		if BaseValue > 3 {
			// Reset state on invalid character
			ForwardHash, ReverseHash = 0, 0
			SmerForwardHash, SmerReverseHash = 0, 0
			SmerPosition = 0
			ValidBases = 0
			continue
		}

		// Rolling forward and reverse K-mer hashes
		ForwardHash = ((ForwardHash << 2) | BaseValue) & ((1 << (2 * K)) - 1)
		Comp := uint64(BaseToInt[Complement[BaseChar]])
		ReverseHash = (ReverseHash >> 2) | (Comp << (2 * (K - 1)))
		ValidBases++

		// Rolling S-mer hashes
		SmerForwardHash = ((SmerForwardHash << 2) | BaseValue) & ((1 << (2 * S)) - 1)
		SmerReverseHash = (SmerReverseHash >> 2) | (Comp << (2 * (S - 1)))
		CanonicalSmer := SmerForwardHash
		if SmerReverseHash < SmerForwardHash {
			CanonicalSmer = SmerReverseHash
		}

		// Store canonical S-mer in ring buffer
		if ValidBases >= S {
			SmerBuffer[SmerPosition%(K-S+1)] = CanonicalSmer
			SmerPosition++

			// When we have filled a full K-mer window
			if ValidBases >= K {
				windowStart := SmerPosition - (K - S + 1)
				minVal := SmerBuffer[0]
				minIdx := 0
				for j := 1; j < K-S+1; j++ {
					idx := (windowStart + j) % (K - S + 1)
					if SmerBuffer[idx] < minVal {
						minVal = SmerBuffer[idx]
						minIdx = j
					}
				}

				if minIdx == T {
					Canonical := ForwardHash
					if ReverseHash < ForwardHash {
						Canonical = ReverseHash
					}
					SyncmerList = append(SyncmerList, SyncmerType{
						Syncmer:  Canonical,
						Position: HeaderID,
					})
				}
			}
		}
	}

	return SyncmerList
}

func SyncmeriseWithWorkerPool(FastaFilePath string, K int, S int, T int, SequenceHashMap map[string]uint64, numWorkers int) []SyncmerType {
	var SyncmerList []SyncmerType
	FastaMap := Fasta.ReadFastaFile(FastaFilePath)

	type Job struct {
		Fasta    Fasta.Fasta
		HeaderID uint64
	}

	jobs := make(chan Job)
	results := make(chan []SyncmerType)

	// Workers
	for w := 0; w < numWorkers; w++ {
		go func() {
			for job := range jobs {
				localSyncmers := SyncmeriseCanonically(job.Fasta, K, S, T, job.HeaderID, nil)
				results <- localSyncmers
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
		localSyncmers := <-results
		SyncmerList = append(SyncmerList, localSyncmers...)
	}
	return SyncmerList
}

func ConvertSyncmerListToKSyncmerMap(SyncmerList []SyncmerType) map[uint64][]uint64 {
	SyncmerMap := make(map[uint64][]uint64)
	for _, Syncmer := range SyncmerList {
		SyncmerMap[Syncmer.Syncmer] = append(SyncmerMap[Syncmer.Syncmer], Syncmer.Position)
	}
	return SyncmerMap
}

func ConvertSyncmerListToKSyncmerMapIgnoringInternalDuplication(SyncmerList []SyncmerType) map[uint64][]uint64 {
	// Ignore marker duplication inside a sequence
	SyncmerMap := make(map[uint64][]uint64)
	for _, Syncmer := range SyncmerList {
		if !slices.Contains(SyncmerMap[Syncmer.Syncmer], Syncmer.Position) {
			SyncmerMap[Syncmer.Syncmer] = append(SyncmerMap[Syncmer.Syncmer], Syncmer.Position)
		}
	}
	return SyncmerMap
}

func FilterSyncmerMap(SyncmerMap map[uint64][]uint64, MinMCount int, MaxMCount int) map[uint64][]uint64 {
	FilteredSyncmerMap := make(map[uint64][]uint64)
	for Syncmer, Positions := range SyncmerMap {
		Count := len(Positions)
		if Count >= MinMCount && Count <= MaxMCount {
			FilteredSyncmerMap[Syncmer] = Positions
		}
	}
	return FilteredSyncmerMap
}

func GenerateSequenceMap(SyncmerMap map[uint64][]uint64) map[uint64]Sequence.Sequence {
	SequenceMap := make(map[uint64]Sequence.Sequence)
	for Syncmer, SequenceList := range SyncmerMap {
		for _, Sequence := range SequenceList {
			info := SequenceMap[Sequence]
			info.Markers = append(info.Markers, Syncmer)
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

func SampleSyncmers(SyncmerMap map[uint64][]uint64, DesiredNumberToRetain int) map[uint64][]uint64 {
	if DesiredNumberToRetain > len(SyncmerMap) {
		DesiredNumberToRetain = len(SyncmerMap)
	}

	// Extract keys
	keys := make([]uint64, 0, len(SyncmerMap))
	for k := range SyncmerMap {
		keys = append(keys, k)
	}

	// Shuffle keys
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	r.Shuffle(len(keys), func(i, j int) {
		keys[i], keys[j] = keys[j], keys[i]
	})

	// Take the first N shuffled keys
	SampledSyncmerMap := make(map[uint64][]uint64)
	for _, k := range keys[:DesiredNumberToRetain] {
		SampledSyncmerMap[k] = SyncmerMap[k]
	}
	return SampledSyncmerMap
}

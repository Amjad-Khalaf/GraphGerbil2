package Kmers

import (
	"GraphGerbil2/Fasta"
	"GraphGerbil2/Sequence"
)

type KmerType struct {
	Kmer     uint64
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

func KmeriseCanonically(Sequence Fasta.Fasta, FilteredKmerList map[uint64]struct{}, K int, HeaderID uint64, KmerList []KmerType) []KmerType {
	// Maintain a hash for forward Kmer, and one for ReverseComplement Kmer
	// Check if the canonical kmer is found in Filtered Tabex List
	// If so, store Canonical Kmer in Kmers slice as KmerType

	SequenceLength := len(Sequence.Sequence)
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
			if _, exists := FilteredKmerList[Canonical]; exists {
				KmerList = append(KmerList, KmerType{Kmer: Canonical, Position: HeaderID})
			}
		}
	}
	return KmerList
}

func KmeriseCanonicallyAllFastaFile(FastaFilePath string, K int, FilteredKmerList map[uint64]struct{}, SequenceHashMap map[string]uint64) []KmerType {
	var KmerList []KmerType
	FastaMap := Fasta.ReadFastaFile(FastaFilePath)
	for _, Fasta := range FastaMap {
		HashedSequenceName := SequenceHashMap[Fasta.Header]
		KmerList = KmeriseCanonically(Fasta, FilteredKmerList, K, HashedSequenceName, KmerList)
	}
	return KmerList
}

func KmeriseWithWorkerPool(FastaFilePath string, K int, FilteredKmerList map[uint64]struct{}, SequenceHashMap map[string]uint64, numWorkers int) []KmerType {
	var KmerList []KmerType
	FastaMap := Fasta.ReadFastaFile(FastaFilePath)

	type Job struct {
		Fasta    Fasta.Fasta
		HeaderID uint64
	}

	jobs := make(chan Job)
	results := make(chan []KmerType)

	// Workers
	for w := 0; w < numWorkers; w++ {
		go func() {
			for job := range jobs {
				localKmers := KmeriseCanonically(job.Fasta, FilteredKmerList, K, job.HeaderID, nil)
				results <- localKmers
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
		localKmers := <-results
		KmerList = append(KmerList, localKmers...)
	}
	return KmerList
}

func ConvertKmerListToKmerMap(KmerList []KmerType) map[uint64][]uint64 {
	KmerMap := make(map[uint64][]uint64)
	for _, Kmer := range KmerList {
		KmerMap[Kmer.Kmer] = append(KmerMap[Kmer.Kmer], Kmer.Position)
	}
	return KmerMap
}

func GenerateSequenceMap(KmerMap map[uint64][]uint64) map[uint64]Sequence.Sequence {
	SequenceMap := make(map[uint64]Sequence.Sequence)
	for Kmer, SequenceList := range KmerMap {
		for _, Sequence := range SequenceList {
			info := SequenceMap[Sequence]
			info.Markers = append(info.Markers, Kmer)
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

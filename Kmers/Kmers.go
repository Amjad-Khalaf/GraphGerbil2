package Kmers

import (
	"GraphGerbil/Fasta"
	"GraphGerbil/Sequence"
	"sort"
	"strconv"
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

func KmeriseCanonically(Sequence Fasta.Fasta, K int, HeaderID uint64, KmerList []KmerType) []KmerType {
	// Maintain a hash for forward Kmer, and one for ReverseComplement Kmer
	// Store Canonical Kmer in Kmers slice as KmerType
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
			KmerList = append(KmerList, KmerType{Kmer: Canonical, Position: HeaderID})
		}
	}
	return KmerList
}

func KmeriseCanonicallyAllFastaFile(FastaFilePath string, K int) ([]KmerType, map[uint64]string) {
	var KmerList []KmerType
	FastaMap := Fasta.ReadFastaFile(FastaFilePath)
	HeaderToIDMap := make(map[uint64]string)
	SequenceNumber := uint64(0)
	for _, Fasta := range FastaMap {
		HeaderToIDMap[SequenceNumber] = Fasta.Header
		KmerList = KmeriseCanonically(Fasta, K, SequenceNumber, KmerList)
	}
	return KmerList, HeaderToIDMap
}

// Update sorting functions so KmerList []KmerType can be sorted
type KmerSlice []KmerType

func (KmerList KmerSlice) Len() int {
	return len(KmerList) // how many KmerType items are in the slice
}
func (KmerList KmerSlice) Less(i, j int) bool {
	return KmerList[i].Kmer < KmerList[j].Kmer // sort by the Kmer field (a uint64)
}
func (KmerList KmerSlice) Swap(i, j int) {
	KmerList[i], KmerList[j] = KmerList[j], KmerList[i] // swap the positions of two elements
}

func FilterKmerList(KmerList []KmerType, MinimumCountInclusive int, MaximumCountInclusive int) []KmerType {
	sort.Sort(KmerSlice(KmerList))
	var FilteredKmerList []KmerType
	for i := 0; i < len(KmerList); {
		current := KmerList[i]
		j := i + 1
		for j < len(KmerList) && KmerList[j].Kmer == current.Kmer {
			j++
		}
		count := j - i
		if count >= MinimumCountInclusive && count <= MaximumCountInclusive {
			FilteredKmerList = append(FilteredKmerList, KmerList[i:j]...)
		}
		i = j
	}
	return FilteredKmerList
}

func ConvertKmerListToKmerMap(KmerList []KmerType, HeaderToIDMap map[uint64]string) map[string][]string {
	KmerMap := make(map[string][]string)
	for _, Kmer := range KmerList {
		KmerStr := strconv.FormatUint(Kmer.Kmer, 10) // convert uint64 to base-10 string
		KmerMap[KmerStr] = append(KmerMap[KmerStr], HeaderToIDMap[Kmer.Position])
	}
	return KmerMap
}

func GenerateSequenceMap(KmerMap map[string][]string) map[string]Sequence.Sequence {
	SequenceMap := make(map[string]Sequence.Sequence)
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

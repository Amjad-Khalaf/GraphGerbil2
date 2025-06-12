package FASTK

import (
	"bufio"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

var BaseToInt = [256]uint8{
	'A': 0, 'C': 1, 'G': 2, 'T': 3,
	'a': 0, 'c': 1, 'g': 2, 't': 3,
}

var Complement = map[byte]byte{
	'A': 'T',
	'C': 'G',
	'G': 'C',
	'T': 'A',
	'a': 't',
	'c': 'g',
	'g': 'c',
	't': 'a',
}

func ReadTabexFileFiltered(TabexFilePath string, K int, MinKmerCount int, MaxKmerCount int, NumOfKmersToRetrieve int) map[uint64]struct{} {
	NumerOfCapturedKmers := 0
	FilteredKmerList := make(map[uint64]struct{})
	File, err := os.Open(TabexFilePath)
	if err != nil {
		log.Fatal(err)
	}
	defer File.Close()

	Reader := bufio.NewReader(File)

	for {
		Line, err := Reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		Line = strings.TrimSpace(Line)
		Fields := strings.Split(Line, "\t")
		if len(Fields) < 2 {
			panic("File is not tab-separated, rerun with Tabex -A")
		}

		KmerString := Fields[0]
		count, err := strconv.Atoi(Fields[1])
		if err != nil {
			log.Fatalf("Invalid count value: %v", err)
		}

		if count < MinKmerCount || count > MaxKmerCount {
			continue // Skip kmers outside the desired range
		}

		// Canonical hash computation
		var ForwardHash, ReverseHash uint64 = 0, 0
		Valid := true
		for i := 0; i < K; i++ {
			if i >= len(KmerString) {
				Valid = false
				break
			}
			Base := KmerString[i]
			BaseValue := BaseToInt[Base]
			if BaseValue > 3 {
				Valid = false
				break
			}
			ForwardHash = ((ForwardHash << 2) | uint64(BaseValue)) & ((1 << (2 * K)) - 1)

			CompBase := Complement[Base]
			CompVal := BaseToInt[CompBase]
			ReverseHash = (ReverseHash >> 2) | (uint64(CompVal) << (2 * (K - 1)))
		} // Move Forward Hash and Reverse Hash Masks outside this loop

		if !Valid {
			continue // skip invalid k-mers
		}

		Canonical := ForwardHash
		if ReverseHash < ForwardHash {
			Canonical = ReverseHash
		}
		FilteredKmerList[Canonical] = struct{}{}
		NumerOfCapturedKmers++
		if NumerOfCapturedKmers == NumOfKmersToRetrieve {
			return FilteredKmerList
		}
	}

	return FilteredKmerList
}

func ReadTabexFile(TabexFilePath string, K int, MinKmerCount int, MaxKmerCount int) map[uint64]struct{} {
	FilteredKmerList := make(map[uint64]struct{})
	File, err := os.Open(TabexFilePath)
	if err != nil {
		log.Fatal(err)
	}
	defer File.Close()

	Reader := bufio.NewReader(File)

	for {
		Line, err := Reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		Line = strings.TrimSpace(Line)
		Fields := strings.Split(Line, "\t")
		if len(Fields) < 2 {
			panic("File is not tab-separated, rerun with Tabex -A")
		}

		KmerString := Fields[0]
		count, err := strconv.Atoi(Fields[1])
		if err != nil {
			log.Fatalf("Invalid count value: %v", err)
		}

		if count < MinKmerCount || count > MaxKmerCount {
			continue // Skip kmers outside the desired range
		}

		// Canonical hash computation
		var ForwardHash, ReverseHash uint64 = 0, 0
		Valid := true
		for i := 0; i < K; i++ {
			if i >= len(KmerString) {
				Valid = false
				break
			}
			Base := KmerString[i]
			BaseValue := BaseToInt[Base]
			if BaseValue > 3 {
				Valid = false
				break
			}
			ForwardHash = ((ForwardHash << 2) | uint64(BaseValue)) & ((1 << (2 * K)) - 1)

			CompBase := Complement[Base]
			CompVal := BaseToInt[CompBase]
			ReverseHash = (ReverseHash >> 2) | (uint64(CompVal) << (2 * (K - 1)))
		}

		if !Valid {
			continue // skip invalid k-mers
		}

		Canonical := ForwardHash
		if ReverseHash < ForwardHash {
			Canonical = ReverseHash
		}
		FilteredKmerList[Canonical] = struct{}{}
	}

	return FilteredKmerList
}

package Fai

import (
	"bufio"
	"io"
	"log"
	"os"
	"slices"
	"strconv"
	"strings"
)

func ReadFaiFile(FaiPath string, SequenceHashMap map[string]uint64) map[uint64]uint64 {
	// Returns:
	// SequenceLengthMap [SequenceName as uint64 from the BUSCO-based Hash Map]Length as uint64

	SequenceLengthMap := make(map[uint64]uint64)

	File, err := os.Open(FaiPath)
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
		i, err := strconv.Atoi(Fields[1])
		if err != nil {
			panic(err)
		}
		SequenceID := SequenceHashMap[Fields[0]]
		SequenceLengthMap[SequenceID] = uint64(i)
	}
	return SequenceLengthMap
}

func SortSequencesBySize(SequenceLengthMap map[uint64]uint64) []uint64 {
	var SequenceList []uint64
	var LengthList []uint64
	for Sequence, Length := range SequenceLengthMap {
		SequenceList = append(SequenceList, Sequence)
		LengthList = append(LengthList, Length)
	}
	var SortedSequenceList []uint64
	SortedLengthList := slices.Clone(LengthList)
	slices.Sort(SortedLengthList)
	slices.Reverse(SortedLengthList)
	var OriginalIndex int
	var UsedIndeces []int
	for _, Length := range SortedLengthList {
		for index, value := range LengthList {
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

func MeasureAssemblySpan(SequenceList []uint64, SequenceLengthMap map[uint64]uint64) uint64 {
	var AssemblySpan uint64 = 0
	for _, Sequence := range SequenceList {
		AssemblySpan = AssemblySpan + SequenceLengthMap[Sequence]
	}
	return AssemblySpan
}

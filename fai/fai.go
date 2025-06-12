package fai

import (
	"bufio"
	"io"
	"log"
	"os"
	"slices"
	"strconv"
	"strings"
)

func ReadFaiFile(FaiPath string, SequenceMap map[string]int) map[string]int {
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
		SequenceMap[Fields[0]] = i
	}
	return SequenceMap
}

func SortSequencesBySize(SequenceMap map[string]int) []string {
	var SequenceList []string
	var LengthList []int
	for Sequence, Length := range SequenceMap {
		SequenceList = append(SequenceList, Sequence)
		LengthList = append(LengthList, Length)
	}
	var SortedSequenceList []string
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

func MeasureAssemblySpan(SequenceList []string, SequenceMap map[string]int) int {
	AssemblySpan := 0
	for _, Sequence := range SequenceList {
		AssemblySpan = AssemblySpan + SequenceMap[Sequence]
	}
	return AssemblySpan
}

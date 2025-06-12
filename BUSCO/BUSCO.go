package BUSCO

import (
	"GraphGerbil/Sequence"
	"bufio"
	"io"
	"log"
	"os"
	"slices"
	"strings"
)

func GenerateSequenceMap(BUSCOFullTable string) map[string]Sequence.Sequence {
	// Returns [Sequence Name] SequenceStruct
	SequenceMap := make(map[string]Sequence.Sequence)
	File, err := os.Open(BUSCOFullTable)
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
		if !strings.Contains(Fields[0], "#") {
			if Fields[1] != "Missing" {
				info := SequenceMap[Fields[2]]
				info.Markers = append(info.Markers, Fields[0])
				info.NumofMarkers++
				SequenceMap[Fields[2]] = info
			}
		}
	}
	return SequenceMap
}

func GenerateBUSCOMap(BUSCOFullTable string) (map[string][]string, int) {
	// Returns [BUSCO] [Sequences which posses said BUSCO]
	BUSCOMap := make(map[string][]string)
	var AllSequenceList []string
	File, err := os.Open(BUSCOFullTable)
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
		if !strings.HasPrefix(Fields[0], "#") {
			if Fields[1] != "Missing" {
				info := BUSCOMap[Fields[0]]
				info = append(info, Fields[2])
				if !slices.Contains(AllSequenceList, Fields[2]) {
					AllSequenceList = append(AllSequenceList, Fields[2])
				}
				BUSCOMap[Fields[0]] = info
			}
		}
	}
	return BUSCOMap, len(AllSequenceList)
}

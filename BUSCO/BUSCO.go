package BUSCO

import (
	"GraphGerbil2/Sequence"
	"bufio"
	"io"
	"log"
	"os"
	"strings"
)

func ReadBUSCOFullTable(BUSCOFullTable string) (map[uint64]Sequence.Sequence, map[string]uint64, map[string]uint64, map[uint64][]uint64) {
	// Returns:
	// 1. SequenceMap [SequenceName as uint64]SequenceStruct
	// 2. SequenceHashMap [SequenceName as string]uint64
	// 3. BUSCOHashMap [BUSCOName as string]uint64
	// 4. BUSCOMap [BUSCOName as uint64][Sequences with BUSCO as uint64]

	// Track sequences for hashing purposes
	var SequenceCount uint64 = 0
	SequenceHashMap := make(map[string]uint64)

	// Track BUSCOs for hashing purposes
	var BUSCOCount uint64 = 0
	BUSCOHashMap := make(map[string]uint64)

	// Sequence Map
	SequenceMap := make(map[uint64]Sequence.Sequence)

	// BUSCO Map
	BUSCOMap := make(map[uint64][]uint64)

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
				SequenceID, ok := SequenceHashMap[Fields[2]]
				if !ok {
					SequenceCount++
					SequenceID = SequenceCount
					SequenceHashMap[Fields[2]] = SequenceID
				}
				info := SequenceMap[SequenceID]
				BUSCOID, ok := BUSCOHashMap[Fields[0]]
				if !ok {
					BUSCOCount++
					BUSCOID = BUSCOCount
					BUSCOHashMap[Fields[0]] = BUSCOID
				}

				// Add information to Sequence Map
				info.Markers = append(info.Markers, BUSCOID)
				info.NumofMarkers++
				SequenceMap[SequenceID] = info

				// Add information to BUSCO Map
				BUSCOMap[BUSCOID] = append(BUSCOMap[BUSCOID], SequenceID)
			}
		}
	}
	return SequenceMap, SequenceHashMap, BUSCOHashMap, BUSCOMap
}

package Fasta

import (
	"bufio"
	"io"
	"log"
	"os"
	"strings"
)

type Fasta struct {
	Header   string
	Sequence string
}

type BufferFasta struct {
	Header   string
	Sequence strings.Builder
}

func ReadFastaFile(path string) map[string]Fasta {
	file, err := os.Open(path)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	var header string
	BufferRecords := make(map[string]BufferFasta)
	Reader := bufio.NewReader(file)
	for {
		line, err := Reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		line = strings.TrimSpace(line)
		if line[0] == '>' {
			header = line[1:]
			BufferRecords[header] = BufferFasta{Header: header}
		} else {
			bufferrecord := BufferRecords[header]
			bufferrecord.Sequence.WriteString(line)
			BufferRecords[header] = bufferrecord
		}
	}
	FastaMap := make(map[string]Fasta)
	for key, value := range BufferRecords {
		FastaMap[key] = Fasta{
			Header:   value.Header,
			Sequence: value.Sequence.String(),
		}
	}
	return FastaMap
}

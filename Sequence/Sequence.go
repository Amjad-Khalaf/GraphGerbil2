package Sequence

import (
	"slices"
)

type Sequence struct {
	Markers      []uint64
	NumofMarkers int
}

type Score struct {
	// Tracks single copy and duplication completeness %s of Markers
	RelativeSingleCopyPercentage float32
	RelativeDuplicatePercentage  float32
	ActualSingleCopyPercentage   float32
	ActualDuplicatePercentage    float32
	DuplicateMarkerNum           int
	SingleCopyMarkerNum          int
	TotalMarkers                 int
	MissingMarkers               int
}

func ScoreMarkers(MarkerList []uint64, TotalMarkerNum int) Score {
	DuplicateNum := 0
	SingleCopyNum := 0
	ScreenedMarkers := make(map[uint64]int)
	for _, Marker := range MarkerList {
		ScreenedMarkers[Marker]++
	}
	for _, count := range ScreenedMarkers {
		if count > 1 {
			DuplicateNum++
		} else if count == 1 {
			SingleCopyNum++
		}
	}
	RelativeTotalMarkerNum := DuplicateNum + SingleCopyNum
	var MarkerScore Score
	MarkerScore.RelativeSingleCopyPercentage = (float32(SingleCopyNum) / float32(RelativeTotalMarkerNum)) * 100
	MarkerScore.RelativeDuplicatePercentage = (float32(DuplicateNum) / float32(RelativeTotalMarkerNum)) * 100
	MarkerScore.ActualSingleCopyPercentage = (float32(SingleCopyNum) / float32(TotalMarkerNum)) * 100
	MarkerScore.ActualDuplicatePercentage = (float32(DuplicateNum) / float32(TotalMarkerNum)) * 100
	MarkerScore.SingleCopyMarkerNum = SingleCopyNum
	MarkerScore.DuplicateMarkerNum = DuplicateNum
	return MarkerScore
}

func JoinTwoSequence(Sequence1Markers []uint64, Sequence2Markers []uint64, TotalMarkerNum int) ([]uint64, Score) {
	JoinedMarkerList := slices.Concat(Sequence1Markers, Sequence2Markers)
	ComparisonScore := ScoreMarkers(JoinedMarkerList, TotalMarkerNum)
	return JoinedMarkerList, ComparisonScore
}

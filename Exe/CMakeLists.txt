SET( CONFOREG_EXE_ITK_SRC
#  seed
#  MedianFiltering
  ResampleImage
  BinaryThresholdFiltering
  DistanceMapCalculation
  KeepRangeFiltering
  DilateDistanceMap
  LocalMaximaExtraction
  SeedsFilteringNuclei
  SeedsGrouping
  )

FOREACH( var ${CONFOREG_EXE_ITK_SRC} )
  ADD_EXECUTABLE( ${var} ${var} )
  TARGET_LINK_LIBRARIES( ${var}
  ITKIO
  ITKCommon
  ITKBasicFilters
  ITKAlgorithms
  ITKStatistics
  ITKNumerics
  )
ENDFOREACH( var ${CONFOREG_EXE_ITK_SRC} )



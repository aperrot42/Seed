/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianSmoothed3DToMembranenessMeasureImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/04/01 21:19:46 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkImageRegionConstIterator.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkListSample.h"
#include <iostream>
#include <queue>
#include <algorithm>

template < class T1, class T2 >
class compareFirst
{
public:
  inline bool operator() (const std::pair<T1,T2>& l, const std::pair<T1,T2>& r)
   {
     return l.first < r.first;
   }
};

template < class T1, class T2 >
class equalSecond
{
public:
  bool operator() (const std::pair<T1,T2>& l, const std::pair<T1,T2>& r)
   {
     return l.second == r.second;
   }
};


int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "seedsImage(itkimage) outputImage(itkimage) " << std::endl
              << "smallestNucleiRadius(double)" << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const unsigned int Dimension = 3;
  // Declare the types of the images
  typedef float       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;
  typedef float       OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension>   OutputImageType;

  // input reader
  typedef itk::ImageFileReader< InputImageType  > ReaderType;

  // Output writer
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  typedef itk::ImageRegionConstIterator< InputImageType > IteratorType;


  typedef itk::Vector< InputImageType::PointValueType, Dimension > MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;

  typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef TreeGeneratorType::KdTreeType TreeType;
  typedef TreeType::NearestNeighbors NeighborsType;
  typedef TreeType::KdTreeNodeType NodeType;
  typedef itk::Statistics::EuclideanDistanceMetric< MeasurementVectorType > DistanceMetricType;

  typedef std::pair< InputImageType::PixelType, InputImageType::PointType > pointPairType;
  typedef std::list< pointPairType > pointListQueueType;


  float   m_smallestRadius = atof(argv[3]);

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  //reader->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );


  reader->GetOutput()->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
  reader->Update();

  pointPairType pixelValuePair;
  pointListQueueType pointsQueue;

  IteratorType inputIterator( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );

  std::cout << "Extracting seeds from image" << std::endl;
  SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize( Dimension );
  MeasurementVectorType mv;

  InputImageType::IndexType pixelIndex;
  InputImageType::PointType pointPosition;
  InputImageType::SpacingType inputSpacing = reader->GetOutput()->GetSpacing();
  InputImageType::PointType inputOrigin = reader->GetOutput()->GetOrigin();
  for ( inputIterator.GoToBegin(); !inputIterator.IsAtEnd(); ++inputIterator)
  {
    if (inputIterator.Get())
      {
      // generate measurement vector from index
      pixelIndex = inputIterator.GetIndex();
      pointPosition[0] = pixelIndex[0]*inputSpacing[0]-inputOrigin[0];
      pointPosition[1] = pixelIndex[1]*inputSpacing[1]-inputOrigin[1];
      pointPosition[2] = pixelIndex[2]*inputSpacing[2]-inputOrigin[2];

      mv[0] = pointPosition[0];
      mv[1] = pointPosition[1];
      mv[2] = pointPosition[2];

      // push to itk sample for generating kdtree
      sample->PushBack( mv );
      // push to shorted list
      pixelValuePair.first = inputIterator.Get();
      pixelValuePair.second =  pointPosition;
      std::cout << pixelValuePair.first << std::endl;
      pointsQueue.push_front(pixelValuePair);
      }
    }

  pointsQueue.sort( compareFirst
                    < InputImageType::PixelType,
                      InputImageType::PointType >() );
 pointsQueue.reverse();

  std::cout << "Generating KdTree" << std::endl;
  TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample( sample );
  treeGenerator->SetBucketSize( 16 );
  treeGenerator->Update();
  TreeType::Pointer tree = treeGenerator->GetOutput();
  NodeType* root = tree->GetRoot();
  if ( root->IsTerminal() )
    {
    std::cout << "Root node is a terminal node." << std::endl;
    }
  else
    {
    std::cout << "Root node is not a terminal node." << std::endl;
    }

  TreeType::InstanceIdentifierVectorType neighbors;


  pointListQueueType::iterator pointListIterator;
  pointListQueueType::iterator pointListIterator2;
  MeasurementVectorType    pointToDeleteMeasurementVector;
  DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
  DistanceMetricType::OriginType origin( Dimension );
  MeasurementVectorType queryPoint;

  pointListQueueType finalPoints;

  std::cout << "Filtering seeds" << std::endl;

  unsigned int pointsToVisit = pointsQueue.size();
  pointListIterator = pointsQueue.begin();
  while ( (pointListIterator != pointsQueue.end()) && pointsToVisit != 0 )
    {
    // create the query point
    pointPosition = pointListIterator->second;
    for ( unsigned int i = 0 ; i < Dimension ; ++i )
      {
      queryPoint[i] = pointPosition[i];
      origin[i] = queryPoint[i];
      }
    //start the query
    tree->Search( queryPoint, m_smallestRadius, neighbors );

    std::cout << "point priority : " << pointListIterator->first
              << " points remaining : " << pointsToVisit
              << std::endl;
    /*
    std::cout << "kd-tree knn search result:" << std::endl
              << "query point = " << queryPoint << std::endl
              << "distance max = " << m_smallestRadius << std::endl;

    distanceMetric->SetOrigin( origin );
    std::cout << "  measurement vector : distance" << std::endl;
    */
    // for all query results
    for ( unsigned int i = 0 ; i < neighbors.size() ; ++i )
      {
      /*
      // for all query results, display distances
      pointToDeleteMeasurementVector = tree->GetMeasurementVector(
          neighbors[i] );
      std::cout << "  "
                << tree->GetMeasurementVector( neighbors[i] )
                << " : "
                << distanceMetric->Evaluate( pointToDeleteMeasurementVector ) << std::endl;
      */
      //delete seeds corresponding to description
      for ( pointListIterator2 = pointsQueue.begin() ;
            pointListIterator2 != pointsQueue.end(); ++pointListIterator2 )
        {
        if (  (pointListIterator2 != pointListIterator)
           && (pointListIterator2->second[0]==(pointToDeleteMeasurementVector[0]) )
           && (pointListIterator2->second[1]==(pointToDeleteMeasurementVector[1]) )
           && (pointListIterator2->second[2]==(pointToDeleteMeasurementVector[2]) ) )
          {
          //std::cout << "    erase : " << pointListIterator2->second << std::endl;
          pointsQueue.erase(pointListIterator2);
          --pointsToVisit;
          }
        }

      }

    // store current query point in the final points vector
    finalPoints.push_front( pointPairType(pointListIterator->first,pointListIterator->second) );
    --pointsToVisit;
    // delete current query point from the point query list
    pointListIterator = pointsQueue.erase(pointListIterator);
    std::cout << std::endl;
    }



  pointListQueueType::iterator finalPointsIterator;


  OutputImageType::IndexType outputPixelIndex;
  OutputImageType::PointType outputPointPosition;
  OutputImageType::SpacingType outputSpacing = reader->GetOutput()->GetSpacing();
  OutputImageType::PointType outputOrigin = reader->GetOutput()->GetOrigin();

  OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer(itk::NumericTraits< OutputPixelType >::Zero);



  std::cout << "writing the final " << finalPoints.size()
            <<" points to an image" << std::endl;

  for ( finalPointsIterator = finalPoints.begin();
        finalPointsIterator != finalPoints.end();
        ++finalPointsIterator )
    {
    outputPointPosition = finalPointsIterator->second;

    for (unsigned int i = 0; i<Dimension;++i)
      {
      outputPixelIndex[i] = (unsigned int) outputPointPosition[i]/inputSpacing[i];
      }
    outputImage->SetPixel(outputPixelIndex,finalPointsIterator->first);
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput ( outputImage );

  std::cout << "writing output image" << std::endl;
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}


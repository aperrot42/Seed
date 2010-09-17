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
     return l.first > r.first;
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

  typedef itk::ImageRegionConstIterator< InputImageType > IteratorType;


  typedef itk::Vector< InputImageType::PointValueType, Dimension > MeasurementVectorType;

  typedef std::pair< InputImageType::PixelType, InputImageType::PointType > pointPairType;
  typedef std::list< pointPairType > pointListQueueType;


  float   m_smallestRadius = atof(argv[3]);

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  //reader->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );


  reader->GetOutput()->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
  reader->Update();

  InputImageType::Pointer inputImage = reader->GetOutput();

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
      // TransformToPhysicalPoint
      inputImage->TransformIndexToPhysicalPoint(pixelIndex,pointPosition);
      /*
      pointPosition[0] = pixelIndex[0]*inputSpacing[0]-inputOrigin[0];
      pointPosition[1] = pixelIndex[1]*inputSpacing[1]-inputOrigin[1];
      pointPosition[2] = pixelIndex[2]*inputSpacing[2]-inputOrigin[2];
      */
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


  pointListQueueType::iterator seedit;
  std::map< float, InputImagePointType > seeds = seedFilter->seeds;
  std::map< float, InputImagePointType >::iterator loc;
  std::cout << "Write out the seed file" << std::endl;
  //Format :
  //Xpos Ypos Zpos Confidence

  std::fstream outfile;
  outfile.open ( argv[3],std::ios::out );

  FeatureImageType::PointType pt;
  
  seedit = pointsQueue.begin();
  while ( seedit != pointsQueue.end() )
    {
    float val = seedit->first;
    pt = seedit->second;
    for ( unsigned int i = 0; i < Dimension; i++ )
      {
      outfile << pt[i] << ' '; 
      }
    outfile << val;
    outfile << std::endl;
    ++seedit;
    }
  outfile.close();

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


int main(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "pointsImage(itkimage) outputSeedTextfile(txt) " << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const unsigned int Dimension = 3;
  // Declare the types of the images
  typedef float       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;

  // input reader
  typedef itk::ImageFileReader< InputImageType  > ReaderType;

  typedef itk::ImageRegionConstIterator< InputImageType > IteratorType;


  typedef std::pair< InputImageType::PixelType, InputImageType::PointType > pointPairType;
  typedef std::list< pointPairType > pointListQueueType;

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

      // push to shorted list
      pixelValuePair.first = inputIterator.Get();
      pixelValuePair.second =  pointPosition;
      std::cout << pixelValuePair.first << std::endl;
      pointsQueue.push_front(pixelValuePair);
      }
    }


  pointListQueueType::iterator seedit;
  std::cout << "Write out the seed file" << std::endl;
  //Format :
  //Xpos Ypos Zpos Confidence

  std::fstream outfile;
  outfile.open ( argv[2],std::ios::out );

  InputImageType::PointType pt;

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
}

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkNumericTraits.h"
#include "itkImage.h"

int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) "
              << "OutputImage(itkimage) "
              << "smallestNucleiRadius(double)" << std::endl;
    return EXIT_FAILURE;
    }
  // Define the dimension of the images
  const int Dimension = 3;

  // Declare the types of the images
  typedef float       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;
  typedef float       OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension>   OutputImageType;

  // input reader
  typedef itk::ImageFileReader< InputImageType  > ReaderType;

  // Output writer
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  // median
  typedef  itk::MedianImageFilter< InputImageType,  InputImageType >
    MedianFilterType;

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();

  float   m_smallestRadius = atof(argv[3]);
  MedianFilterType::RadiusType structuringRadius;
  for (int i=0; i< Dimension;++i)
    {
    structuringRadius[i] = static_cast<unsigned int>(
                             m_smallestRadius
                             /(reader->GetOutput()->GetSpacing()[i]) );
    }

   MedianFilterType::Pointer medianFilter =  MedianFilterType::New();
   medianFilter->SetRadius( structuringRadius );

  std::cout << "median filtering on a neighborhood of radius : " 
            << medianFilter->GetRadius()
            << std::endl;
  medianFilter->SetInput(reader->GetOutput());
  medianFilter->SetNumberOfThreads(6);
  medianFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( medianFilter->GetOutput() );
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

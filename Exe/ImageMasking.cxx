#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkNumericTraits.h"

int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "inputImage(itkimage) mask(itkimage) "
              << "OutputImage(itkimage)" << std::endl;
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

  // distance function filter
  typedef itk::MultiplyImageFilter< InputImageType, InputImageType, OutputImageType >
    MultiplyFilterType;

  std::cout << "reading inputImage" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );

  std::cout << "reading mask" << std::endl;
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName ( argv[2] );

  std::cout << "multiply seeds by nuclei mask computing" << std::endl;
  MultiplyFilterType::Pointer multiplyFilter
    = MultiplyFilterType::New();
  multiplyFilter->SetInput1(reader->GetOutput());
  multiplyFilter->SetInput2(reader2->GetOutput());
  multiplyFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput ( multiplyFilter->GetOutput() );

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

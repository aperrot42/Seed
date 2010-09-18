#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkNumericTraits.h"
#include "itkImage.h"

int main(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) "
              << std::endl;
    return EXIT_FAILURE;
    }
  // Define the dimension of the images
  const int Dimension = 3;

  // Declare the types of the images
  typedef float       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;
  typedef unsigned char       OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension>   OutputImageType;

  // input reader
  typedef itk::ImageFileReader< InputImageType  > ReaderType;

  // Output writer
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  // rescale
  typedef itk::RescaleIntensityImageFilter< InputImageType,  InputImageType >
    RescaleFilterType;
    
   typedef itk::CastImageFilter< InputImageType,  OutputImageType >
    CastFilterType;

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();
  

  InputImageType::Pointer inputImage = reader->GetOutput();

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();

  rescaleFilter->SetNumberOfThreads(6);
  rescaleFilter->SetInput(inputImage);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetOutputMinimum(0);

  rescaleFilter->Update();

  CastFilterType::Pointer castFilter = CastFilterType::New();

  castFilter->SetNumberOfThreads(6);
  castFilter->SetInput(rescaleFilter->GetOutput());
  castFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( castFilter->GetOutput() );

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
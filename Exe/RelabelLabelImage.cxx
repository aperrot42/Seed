#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include <stdio.h>

#include "itkRelabelComponentImageFilter.h"

int main(int argc, char * argv[])
{

  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) OutputImage(itkimage)" << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const int Dimension = 3;

  // Declare the types of the images
  typedef unsigned int       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;
  typedef InputPixelType OutputPixelType;
  typedef InputImageType OutputImageType;

  // read the input image
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::RelabelComponentImageFilter< InputImageType, OutputImageType > RelabelFilterType;
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  relabelFilter->SetInput(reader->GetOutput());
  relabelFilter->Update();
  
  std::cout << "Number of objects " << relabelFilter->GetNumberOfObjects() << std::endl;

    typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(relabelFilter->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}

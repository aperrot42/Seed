#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include <stdio.h>
#include <math.h>

int main(int argc, char * argv[])
{

  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) OutputImage(itkImage) nucleiradius" << std::endl;
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

  typedef unsigned int BinaryPixelType;
  typedef itk::Image<BinaryPixelType, Dimension> BinaryImageType;


  // define the object type. Here the ShapeLabelObject type
  // is chosen in order to read some attribute related to the shape
  // of the objects (by opposition to the content of the object, with
  // the StatisticsLabelObejct).
  typedef unsigned long LabelType;
  typedef itk::ShapeLabelObject< LabelType, Dimension > LabelObjectType;
  typedef itk::LabelMap< LabelObjectType > LabelMapType;

  // convert the image in a collection of objects
  typedef itk::LabelImageToShapeLabelMapFilter< BinaryImageType, LabelMapType > ConverterType;
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInput( reader->GetOutput() );
  converter->Update();
  
  
  // create output image
  OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->SetSpacing(reader->GetOutput()->GetSpacing());
  outputImage->SetOrigin(reader->GetOutput()->GetOrigin());
  outputImage->Allocate();
  outputImage->SetBufferedRegion(outputImage->GetLargestPossibleRegion() );
  outputImage->FillBuffer(0);
  
  OutputImageType::PointType outPoint;
  OutputImageType::IndexType outPix;

  BinaryImageType::PointType centroid;
  std::fstream outfile;
  outfile.open ( argv[2],std::ios::out );
  LabelMapType::Pointer labelMap = converter->GetOutput();
  for( unsigned int label=1; label<=labelMap->GetNumberOfLabelObjects(); label++ )
    {
    // we don't need a SmartPointer of the label object here, because the reference is kept in
    // in the label map.
    const LabelObjectType * labelObject = labelMap->GetLabelObject( label );
    std::cout << label << "\t" <<  labelObject->GetCentroid();
    centroid = labelObject->GetCentroid();
    outPoint = centroid;
    outputImage->TransformPhysicalPointToIndex(outPoint, outPix);
    outputImage->SetPixel(outPix, label);
    std::cout << centroid[0] << ' ' << centroid[1] << ' '<< centroid[2] << ' ';
std::cout << std::endl;
    }
  outfile.close();

  float nucleiRadius = (atof(argv[3]));
  typedef itk::BinaryBallStructuringElement< OutputPixelType, Dimension >
    BinaryStructuringType;
  typedef itk::GrayscaleDilateImageFilter<OutputImageType,OutputImageType,
    BinaryStructuringType> DilateFilterType;

  DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
  
  BinaryStructuringType  structuringElement;
  BinaryStructuringType::RadiusType ballRadius;
  ballRadius[0] = (unsigned int) (nucleiRadius/(outputImage->GetSpacing()[0]));
  ballRadius[1] = (unsigned int) (nucleiRadius/(outputImage->GetSpacing()[1]));
  ballRadius[2] = (unsigned int) (nucleiRadius/(outputImage->GetSpacing()[2]));
  std::cout << ballRadius << std::endl;
  structuringElement.SetRadius( ballRadius );
  structuringElement.CreateStructuringElement();
  
  dilateFilter->SetKernel( structuringElement );

  dilateFilter->SetInput(outputImage);
  
  dilateFilter->SetNumberOfThreads(8);
  dilateFilter->Update();


  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(dilateFilter->GetOutput());
  writer->Update();
  

  return EXIT_SUCCESS;
}

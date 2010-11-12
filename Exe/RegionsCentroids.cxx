#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkGeometryUtilities.h"
#include "itkBinaryThresholdImageFilter.h"
#include <stdio.h>






int main(int argc, char * argv[])
{

  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) OutputPoints(txt) [largestCellRadius]" << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  const int Dimension = 3;

  // Declare the types of the images
  typedef float       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;

  // read the input image
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef unsigned int BinaryPixelType;
  typedef itk::Image<BinaryPixelType, Dimension> BinaryImageType;
  typedef itk::BinaryThresholdImageFilter< InputImageType, BinaryImageType > BinaryFilterType;
  BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();
  binaryFilter->SetInput(reader->GetOutput());
  binaryFilter->SetInsideValue(1);
  binaryFilter->SetOutsideValue(0);
  binaryFilter->SetLowerThreshold(0.01);
  binaryFilter->SetUpperThreshold(1000);
  binaryFilter->Update();
  

  typedef itk::ConnectedComponentImageFilter< BinaryImageType, BinaryImageType > ConnectedComponentType;
  ConnectedComponentType::Pointer connectedFilter = ConnectedComponentType::New();
  connectedFilter->SetInput(binaryFilter->GetOutput());
  connectedFilter->Update();
  connectedFilter->SetBackgroundValue(0);
  
    typedef itk::ImageFileWriter< BinaryImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "toto.mha" );
  writer->SetInput(connectedFilter->GetOutput());
  writer->Update();
  
  // define the object type. Here the ShapeLabelObject type
  // is chosen in order to read some attribute related to the shape
  // of the objects (by opposition to the content of the object, with
  // the StatisticsLabelObejct).
  typedef unsigned int LabelType;
  typedef itk::ShapeLabelObject< LabelType, Dimension > LabelObjectType;
  typedef itk::LabelMap< LabelObjectType > LabelMapType;

  // convert the image in a collection of objects
  typedef itk::LabelImageToShapeLabelMapFilter< BinaryImageType, LabelMapType > ConverterType;

  ConverterType::Pointer converter = ConverterType::New();
/*
  converter->SetInput( connectedFilter->GetOutput() );
  //converter->SetInputForegroundValue(2);
  //converter->FullyConnectedOn ();
  // update the shape filter, so its output will be up to date
  //converter->SetOutputBackgroundValue (0);


  converter->Update();
  
  
  float LargestCellRadius = 1.;

  if (argc > 3)
    {
    LargestCellRadius = (atof(argv[3]));
    LargestCellRadius *= LargestCellRadius;
    LargestCellRadius *= LargestCellRadius;
    LargestCellRadius *= 3.14*4./3.;
    }
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
    if ((labelObject->GetNumberOfPixels() > 4) || (labelObject->GetPhysicalSize() > LargestCellRadius*5 ))
      {
      outfile << centroid[0] << ' ' << centroid[1] << ' '<< centroid[2] << ' ';

      if (argc > 3)
        {
        std::cout << "\t" << labelObject->GetPhysicalSize()/(LargestCellRadius) << std::endl;
        // confidence decrease if inside of cell too big
        outfile << (labelObject->GetPhysicalSize())/(LargestCellRadius);
        }
      else
        {
        std::cout << "\t" << LargestCellRadius << std::endl;
        // if we have no info on info
        outfile <<LargestCellRadius;
        }
      outfile << std::endl;
      }
    else
      {
      if (labelObject->GetNumberOfPixels() > 4)
        std::cout << "Refused : too small" << std::endl;
      else
        std::cout << "Refused : too big" << std::endl;
      }
    }
  outfile.close();


*/

  return EXIT_SUCCESS;
}

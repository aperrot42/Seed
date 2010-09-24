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
#include "itkNumericTraits.h"
#include "itkDoubleThresholdImageFilter.h"

int main(int argc, char* argv [] )
{
  if ( argc < 7 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "inputImage outputImage" << std::endl
              << "Threshold1 Threshold2 Threshold3 Threshold4" << std::endl;
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

  // binary threshold filter
  typedef  itk::DoubleThresholdImageFilter< InputImageType, InputImageType >
    DoubleThresholdFilterType;

  float Threshold1 = atof(argv[3]);
  float Threshold2 = atof(argv[4]);
  float Threshold3 = atof(argv[5]);
  float Threshold4 = atof(argv[6]);
  
  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();


  std::cout << "Threshold computing" << std::endl;
  std::cout << "Threshold1 =" << Threshold1 << std::endl;
  std::cout << "Threshold2 =" << Threshold2 << std::endl;
  std::cout << "Threshold2 =" << Threshold3 << std::endl;
  std::cout << "Threshold2 =" << Threshold4 << std::endl;

  DoubleThresholdFilterType::Pointer binaryFilter
    = DoubleThresholdFilterType::New();
  binaryFilter->SetInput(reader->GetOutput());
  binaryFilter->SetThreshold1(Threshold1);
  binaryFilter->SetThreshold2(Threshold2);
  binaryFilter->SetThreshold3(Threshold3);
  binaryFilter->SetThreshold4(Threshold4);

  binaryFilter->SetInsideValue (itk::NumericTraits< OutputPixelType >::One );
  binaryFilter->SetOutsideValue (itk::NumericTraits< OutputPixelType >::Zero );
  binaryFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( binaryFilter->GetOutput() );

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


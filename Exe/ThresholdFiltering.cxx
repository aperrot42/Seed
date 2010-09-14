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

int main(int argc, char* argv [] )
{
  if ( argc < 5 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "inputImage outputImage" << std::endl
              << "LowThreshold HighThreshold" << std::endl;
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
  typedef itk::BinaryThresholdImageFilter< InputImageType, InputImageType >
    BinaryFilterType;

  float m_lowThreshold;
  float m_highThreshold;

  m_lowThreshold = atof(argv[3]);
  m_highThreshold = atof(argv[4]);

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();


  std::cout << "Threshold computing" << std::endl;
  std::cout << "LowThreshold =" << m_lowThreshold << std::endl;
  std::cout << "HighThreshold =" << m_highThreshold << std::endl;

  BinaryFilterType::Pointer binaryFilter
    = BinaryFilterType::New();
  binaryFilter->SetInput(reader->GetOutput());
  binaryFilter->SetUpperThreshold(m_highThreshold);
  binaryFilter->SetLowerThreshold(m_lowThreshold);
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


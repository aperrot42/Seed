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
#include "itkThresholdImageFilter.h"
#include "itkNumericTraits.h"

int main(int argc, char* argv [] )
{
  if ( argc < 5 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "inputImage(itkimage) outputImage(itkimage)" << std::endl
              << "LowThreshold(double) HighThreshold(double) "
              << "OutOfRange(double)" << std::endl;
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
  typedef itk::ThresholdImageFilter< InputImageType >
    thresholdFilterType;

  float m_lowThreshold = atof(argv[3]);
  float m_highThreshold = atof(argv[4]);
  float m_outOfRange = atof(argv[5]);

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();

  std::cout << "Threshold computing" << std::endl;
  std::cout << "LowThreshold = " << m_lowThreshold << std::endl;
  std::cout << "HighThreshold = " << m_highThreshold << std::endl;
  std::cout << "OutOfRang value = " << m_outOfRange << std::endl;

  thresholdFilterType::Pointer thresholdFilter
    = thresholdFilterType::New();
  thresholdFilter->SetInput(reader->GetOutput());
  thresholdFilter->ThresholdOutside( m_lowThreshold, m_highThreshold );
  thresholdFilter->SetOutsideValue( m_outOfRange);
  thresholdFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( thresholdFilter->GetOutput() );

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


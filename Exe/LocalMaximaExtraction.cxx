
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkValuedRegionalMaximaImageFilter.h"
#include "itkThresholdImageFilter.h"
#include <stdio.h>


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
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) OutputImage(itkimage)" << std::endl
              << "MaximumMaxima(double)" << std::endl;
    return EXIT_FAILURE;
    }

  float m_MaxDistance = atof(argv[3]);

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

  // localmax function filter
  typedef itk::ValuedRegionalMaximaImageFilter< InputImageType, InputImageType >
    LocalMaximaFilterType;

  // threshold filter
  typedef itk::ThresholdImageFilter< OutputImageType>
    ThresholdFilterType;

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();

  std::cout << "Local maxima computing" << std::endl;
  std::cout << "maximum local maxima=" << m_MaxDistance << std::endl;

  LocalMaximaFilterType::Pointer localMaximaFilter
    =  LocalMaximaFilterType::New();
  localMaximaFilter->SetInput(reader->GetOutput());
  localMaximaFilter->Update();


  ThresholdFilterType::Pointer thresholdFilter
      = ThresholdFilterType::New();
  thresholdFilter->SetInput(localMaximaFilter->GetOutput());
  thresholdFilter->ThresholdOutside(itk::NumericTraits< OutputPixelType >::Zero,
                                    m_MaxDistance);
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

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include  "itkNearestNeighborInterpolateImageFunction.h"

int main( int argc, char * argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputImageFile(itkimage)  outputImageFile(itkimage) "
              << " lowerIntensity(uchar) upperIntensity(uchar) Nearest/Linear/Bspline(-1/0/5)" << std::endl;
    return EXIT_FAILURE;
    }

  const     unsigned int    Dimension = 3;

  typedef   unsigned short  InputPixelType;
  typedef   float           InternalPixelType;

  typedef itk::Image< InputPixelType,    Dimension >   InputImageType;
  typedef itk::Image< InternalPixelType, Dimension >   InternalImageType;



  typedef itk::ImageFileReader< InputImageType  >  ReaderType;

  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << excep << std::endl;
    }

  typedef itk::IntensityWindowingImageFilter<
                                  InputImageType,
                                  InternalImageType >  IntensityFilterType;



  IntensityFilterType::Pointer intensityWindowing = IntensityFilterType::New();

  intensityWindowing->SetWindowMinimum( atoi( argv[3] ) );
  intensityWindowing->SetWindowMaximum( atoi( argv[4] ) );

  intensityWindowing->SetOutputMinimum(   0.0 );
  intensityWindowing->SetOutputMaximum( 255.0 ); // floats but in the range of chars.

  intensityWindowing->SetInput( reader->GetOutput() );


  typedef itk::RecursiveGaussianImageFilter<
                                InternalImageType,
                                InternalImageType > GaussianFilterType;


  GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
  GaussianFilterType::Pointer smootherY = GaussianFilterType::New();

  smootherX->SetInput( intensityWindowing->GetOutput() );
  smootherY->SetInput( smootherX->GetOutput() );

  InputImageType::ConstPointer inputImage = reader->GetOutput();

  const InputImageType::SpacingType& inputSpacing = inputImage->GetSpacing();

  const double isoSpacing = vcl_sqrt( inputSpacing[2] * inputSpacing[0] );

  smootherX->SetSigma( isoSpacing*0.3 );
  smootherY->SetSigma( isoSpacing*0.3 );

  smootherX->SetDirection( 0 );
  smootherY->SetDirection( 1 );

  smootherX->SetNormalizeAcrossScale( true );
  smootherY->SetNormalizeAcrossScale( true );

  typedef   float   OutputPixelType;

  typedef itk::Image< OutputPixelType,   Dimension >   OutputImageType;

  typedef itk::ResampleImageFilter<
                InternalImageType, OutputImageType >  ResampleFilterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();


  typedef itk::IdentityTransform< double, Dimension >  TransformType;

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  resampler->SetTransform( transform );


  typedef itk::BSplineInterpolateImageFunction<
    InternalImageType, double >  SplineInterpolatorType;
  SplineInterpolatorType::Pointer splineInterpolator =
    SplineInterpolatorType::New();

  typedef itk::LinearInterpolateImageFunction<
    InternalImageType, double >  LinearInterpolatorType;
  LinearInterpolatorType::Pointer linearInterpolator =
    LinearInterpolatorType::New();

   typedef itk::NearestNeighborInterpolateImageFunction<
    InternalImageType, double >  NearesInterpolatorType;
  NearesInterpolatorType::Pointer nearestInterpolator =
    NearesInterpolatorType::New();

  if ( (atoi( argv[5] )) > 0) //Bspline 0-5
    {
    splineInterpolator->SetSplineOrder(atoi( argv[5] ));
    resampler->SetInterpolator(splineInterpolator);
    }
  else
    {
    if ( (atoi( argv[5] ) ) == 0)// Linear
      {
      resampler->SetInterpolator( linearInterpolator );
      }
    else // nearest (-1)
      {
      resampler->SetInterpolator( nearestInterpolator );
      }
    }


  resampler->SetDefaultPixelValue( 0 );

  OutputImageType::SpacingType spacing;

  spacing[0] = isoSpacing;
  spacing[1] = isoSpacing;
  spacing[2] = isoSpacing;

  resampler->SetOutputSpacing( spacing );

  resampler->SetOutputOrigin( inputImage->GetOrigin() );
  resampler->SetOutputDirection( inputImage->GetDirection() );

  InputImageType::SizeType   inputSize =
                    inputImage->GetLargestPossibleRegion().GetSize();

  typedef InputImageType::SizeType::SizeValueType SizeValueType;

  const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
  const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;

  const double dz = (inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;

  InputImageType::SizeType   size;

  size[0] = static_cast<SizeValueType>( dx );
  size[1] = static_cast<SizeValueType>( dy );
  size[2] = static_cast<SizeValueType>( dz );

  resampler->SetSize( size );
  resampler->SetInput( smootherY->GetOutput() );

  smootherX->SetNumberOfThreads(6);
  smootherY->SetNumberOfThreads(6);
  //splineInterpolator->SetNumberOfThreads(2);
  resampler->SetNumberOfThreads(6);

  resampler->Update();

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[2] );
  writer->SetInput( resampler->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }


  return EXIT_SUCCESS;
}


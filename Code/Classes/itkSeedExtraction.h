/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1634 $  // Revision of last commit
  Date: $Date: 2010-06-10 17:01:47 -0400 (Thu, 10 Jun 2010) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef __itkSeedExtraction_h
#define __itkSeedExtraction_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageToImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkImageRegion.h"
#include "itkRegion.h"
#include "itkIndex.h"
#include "itkSize.h"
#include <map>
#include <vector>
#include <vnl/vnl_vector.h>

#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
template < class TFeatureImage, class TInputImage, class TSegmentImage >
  class ITK_EXPORT SeedExtraction :
  public ImageToImageFilter< TInputImage, TSegmentImage >
{
  public:
    typedef SeedExtraction             Self;
    typedef ImageToImageFilter< TInputImage,TSegmentImage> Superclass;
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;

    itkStaticConstMacro ( ImageDimension, unsigned int,
                          TFeatureImage::ImageDimension );

    /** Method for creation through object factory */
    itkNewMacro ( Self );

    /** Run-time type information */
    itkTypeMacro ( SeedExtraction, ImageToImageFilter );

    /** Display */
    void PrintSelf ( std::ostream& os, Indent indent ) const;

    typedef TFeatureImage                           FeatureImageType;
    typedef typename FeatureImageType::Pointer      FeatureImagePointer;
    typedef typename FeatureImageType::ConstPointer FeatureImageConstPointer;
    typedef typename FeatureImageType::PixelType    FeatureImagePixelType;
    typedef typename FeatureImageType::RegionType   FeatureImageRegionType;
    typedef typename FeatureImageType::SizeType     FeatureImageSizeType;
    typedef typename FeatureImageSizeType::SizeValueType FeatureImageSizeValueType;
    typedef typename FeatureImageType::SpacingType  FeatureImageSpacingType;
    typedef typename FeatureImageType::IndexType    FeatureImageIndexType;
    typedef typename FeatureImageType::PointType    FeatureImagePointType;

    typedef TInputImage                             ImageType;
    typedef typename ImageType::Pointer             ImagePointer;
    typedef typename ImageType::ConstPointer        ImageConstPointer;
    typedef typename ImageType::PixelType           ImagePixelType;
    typedef typename ImageType::RegionType          ImageRegionType;
    typedef typename ImageType::SizeType            ImageSizeType;
    typedef typename ImageSizeType::SizeValueType   ImageSizeValueType;
    typedef typename ImageType::SpacingType         ImageSpacingType;
    typedef typename ImageType::IndexType           ImageIndexType;
    typedef typename ImageIndexType::IndexValueType ImageIndexValueType;
    typedef typename ImageType::PointType           ImagePointType;

    typedef TSegmentImage                           SegmentImageType;
    typedef typename SegmentImageType::Pointer      SegmentImagePointer;
    typedef typename SegmentImageType::ConstPointer SegmentImageConstPointer;
    typedef typename SegmentImageType::IndexType    SegmentImageIndexType;
    typedef typename SegmentImageType::PixelType    SegmentImagePixelType;

    typedef CastImageFilter< ImageType, ImageType >  InputCastType;
    typedef typename InputCastType::Pointer InputCastPointer;
    typedef ImageRegionIterator< ImageType > IteratorType;
    typedef ImageRegionConstIterator< ImageType > ConstIteratorType;
		typedef ImageRegionConstIteratorWithIndex< ImageType > ConstIndexIteratorType;
    typedef ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;
    typedef ImageRegionIterator< SegmentImageType > SegmentIteratorType;
    typedef ImageRegionIteratorWithIndex< SegmentImageType > SegmentIndexIteratorType;

    typedef DiscreteGaussianImageFilter< ImageType, ImageType >
      GaussianFilterType;
    typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

    typedef MinimumMaximumImageCalculator< ImageType > MinMaxCalculatorType;
    typedef typename MinMaxCalculatorType::Pointer     MinMaxCalculatorPointer;

    typedef ConnectedComponentImageFilter< SegmentImageType, SegmentImageType >
      ConnectedComponentFilterType;
    typedef typename ConnectedComponentFilterType::Pointer
      ConnectedComponentFilterPointer;

    typedef RelabelComponentImageFilter< SegmentImageType, SegmentImageType >
      RelabelComponentFilterType;
    typedef typename RelabelComponentFilterType::Pointer
      RelabelComponentFilterPointer;

  	typedef BinaryBallStructuringElement< ImagePixelType, ImageDimension> SEType;
	  typedef GrayscaleDilateImageFilter< ImageType, ImageType, SEType > DilateFilterType;

    typedef vnl_vector< float > vnlVectorType;
    typedef itk::Vector< float, ImageDimension > MeasurementVectorType;
    typedef Statistics::ListSample< MeasurementVectorType > SampleType;
    typedef typename SampleType::Pointer SamplePointer;

    typedef Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
    typedef typename TreeGeneratorType::Pointer TreeGeneratorPointer;
    typedef typename TreeGeneratorType::KdTreeType TreeType;
    typedef typename TreeType::Pointer TreePointer;
    typedef typename TreeType::InstanceIdentifierVectorType TreeInstanceId;
    typedef itk::Statistics::EuclideanDistanceMetric< MeasurementVectorType >
      DistanceMetricType;
    typedef typename DistanceMetricType::Pointer DistanceMetricPointer;

    itkGetConstMacro ( NumberOfSeeds, unsigned int );
    itkSetMacro ( NumberOfSeeds, unsigned int );
    itkGetConstMacro ( LargestCellRadius, double );
    itkSetMacro ( LargestCellRadius, double );
    itkGetConstMacro ( Threshold, double );
    itkSetMacro ( Threshold, double );

    typedef typename std::map< float, ImagePointType >::iterator SeedIteratorType;
    std::map< float, ImagePointType > seeds;
    std::vector< unsigned long >      sizes;

    void SetForeground ( SegmentImagePointer fg )
    {
      m_fgMap = fg;
      this->Modified();
    }

  protected:
    SeedExtraction();
    ~SeedExtraction() {}

    void ImageSeeds ();
    void GenerateData();
  	void RemoveRegion( ImagePointType pt );
    void NonMaximalSeeds();

    double              m_LargestCellRadius;
    double              m_Threshold;
    unsigned int        m_NumberOfSeeds;
	  ImageSpacingType    m_Spacing;
	  ImageSizeType       m_Size;
    SamplePointer   m_Sample;
    TreePointer m_Tree;

    ImagePointer m_Input;
    SegmentImagePointer m_fgMap;
    SegmentImagePointer m_relabelfgMap;

  private:
    SeedExtraction ( Self& );       // intentionally not implemented
    void operator= ( const Self& ); // intentionally not implemented
};

} /* namespace itk */

#include "itkSeedExtraction.txx"
#endif

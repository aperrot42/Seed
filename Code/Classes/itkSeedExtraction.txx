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

#ifndef __itkSeedExtraction_txx
#define __itkSeedExtraction_txx

#include "itkSeedExtraction.h"

namespace itk
{
//  Software Guide : BeginCodeSnippet
template < class TFeatureImage, class TInputImage, class TSegmentImage >
SeedExtraction< TFeatureImage, TInputImage, TSegmentImage >
::SeedExtraction()
{
  m_LargestCellRadius = 5.0;
  m_Threshold = 0.5;

  this->Superclass::SetNumberOfRequiredInputs ( 1 );
  this->Superclass::SetNumberOfRequiredOutputs ( 1 );

  this->Superclass::SetNthOutput ( 0, TSegmentImage::New() );
}


template < class TFeatureImage, class TInputImage, class TSegmentImage >
void
SeedExtraction< TFeatureImage, TInputImage, TSegmentImage >
::NonMaximalSeeds ()
{
  MinMaxCalculatorPointer minMaxCalculator = MinMaxCalculatorType::New();
  minMaxCalculator->SetImage( m_Input );

  ImageIndexType idx;
  ImagePointType pt;
  ImagePixelType max;

  unsigned int circles = seeds.size();

  // Find maxima
  do
    {
    minMaxCalculator->ComputeMaximum();
    max = minMaxCalculator->GetMaximum();
    idx = minMaxCalculator->GetIndexOfMaximum();
    this->m_Input->TransformIndexToPhysicalPoint( idx, pt );

    if ( max < m_Threshold )
      {
      break;
      }

    while ( seeds.find ( max ) != seeds.end() )
      max += 0.001;

    seeds[max] = pt;

    RemoveRegion( pt );
    ++circles;
    } while(circles < m_NumberOfSeeds);
}


template < class TFeatureImage, class TInputImage, class TSegmentImage >
void
SeedExtraction< TFeatureImage, TInputImage, TSegmentImage >
::RemoveRegion( ImagePointType pt )
{
  ImageRegionType region;
  ImageIndexType start, end, idx;
  ImageSizeType sizeOfROI;
  ImageIndexValueType rad;
  this->GetInput()->TransformPhysicalPointToIndex( pt, idx );

  // Remove a black disc from the hough space domain
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    rad = static_cast<ImageIndexValueType>( 1.7*m_LargestCellRadius/m_Spacing[i] );

    if( idx[i] > static_cast<ImageIndexValueType>( rad ) )
      {
      start[i] = idx[i] - static_cast<ImageIndexValueType>( rad );
      }
    else
      {
      start[i] = 0;
      }

    if( idx[i] + rad < static_cast<ImageIndexValueType>( m_Size[i] - 1 ) )
      {
      end[i] = idx[i] + static_cast<ImageIndexValueType>( rad );
      }
    else
      {
      end[i] = m_Size[i] - 1;
      }

    sizeOfROI[i] = end[i] - start[i] + 1;
    }

  region.SetSize( sizeOfROI );
  region.SetIndex( start );

  // Delete the region values
  IndexIteratorType sIt ( m_Input, region );
  for ( sIt.GoToBegin(); !sIt.IsAtEnd(); ++sIt )
    {
    sIt.Set( 0 );
    }
}


template < class TFeatureImage, class TInputImage, class TSegmentImage >
void
SeedExtraction< TFeatureImage, TInputImage, TSegmentImage >
::ImageSeeds ()
{
  // Specify the dilation radius of the structuring mask. It correlates with the
  // minimum separation between the nuclei.
  ImageSpacingType spacing = this->GetInput()->GetSpacing();
  ImageSizeType radius;
  for ( unsigned int j = 0; j < ImageDimension; j++ )
  {
    radius[j] = static_cast<ImageSizeValueType> (
                  0.4*m_LargestCellRadius/spacing[j] );
  }


  ImagePixelType m = 0;
  ImageRegionType region = m_relabelfgMap->GetLargestPossibleRegion();
  IndexIteratorType sIt ( m_Input, region );
  SegmentIteratorType fIt ( m_relabelfgMap, region );
  for ( sIt.GoToBegin(), fIt.GoToBegin();
        !sIt.IsAtEnd(); ++sIt, ++fIt )
    {
    if ( fIt.Get() > 0 ) 
      {
      if ( m < sIt.Get() )
        {
        m = sIt.Get();
        }
      }
    else
      {
      sIt.Set( 0 );
      }
    }
  m_Threshold *= m;

  std::cout << "  Maxima " << m <<" Threshold " << m_Threshold << std::endl;

  SEType sE;
  sE.SetRadius ( radius );
  sE.CreateStructuringElement();

  // Dilate the input volume
  typename DilateFilterType::Pointer grayscaleDilate = DilateFilterType::New();
  grayscaleDilate->SetKernel ( sE );
  grayscaleDilate->SetInput ( this->GetInput() );
  grayscaleDilate->Update();
  ImagePointer dilateImage = grayscaleDilate->GetOutput();

  std::map< float, ImagePointType > allSeeds;
  MeasurementVectorType mv;
  ImagePointType pt;
  ImageIndexType idx;
  float val = 0;
  bool flag;
  // Iterate over the smooth distance field
  IteratorType dIt ( dilateImage, region );
  for ( sIt.GoToBegin(), dIt.GoToBegin(), fIt.GoToBegin();
        !sIt.IsAtEnd();
       ++sIt, ++dIt, ++fIt )
    {
    // If the pixel after dilation has the same value prior to dilation, it is a
    //local maxima! Make sure it lies in the cell foreground.
    if ( ( sIt.Get() == dIt.Get() ) && ( fIt.Get() > 0 ) )
      {
      flag = sIt.Get() > m_Threshold;
      val = sIt.Get() /m;
      if ( flag )
        {
        while ( allSeeds.find ( val ) != allSeeds.end() )
          val += 0.001;

			  idx = sIt.GetIndex();
        dilateImage->TransformIndexToPhysicalPoint( idx, pt );
        allSeeds[val] = pt;
        }
      }
    }

  SeedIteratorType it = allSeeds.begin();
  m_Sample = SampleType::New();
  m_Sample->SetMeasurementVectorSize( ImageDimension );

  while ( it != allSeeds.end() )
    {
    // Insert into a kd-tree here
    pt = it->second;
    for( unsigned int i = 0 ; i < ImageDimension; i++ )
      {
      mv[i] = static_cast<float>( pt[i] );
      }
      m_Sample->PushBack( mv );
    ++it;
    }

  TreeGeneratorPointer m_TreeGenerator = TreeGeneratorType::New();
  m_TreeGenerator->SetSample( m_Sample );
  m_TreeGenerator->SetBucketSize( 4 ); // this value can be modified
  m_TreeGenerator->Update();
  m_Tree = m_TreeGenerator->GetOutput();

  unsigned int numberOfNeighbors = 5;
  TreeInstanceId neighbors;

  float dist;
  DistanceMetricPointer distanceMetric = DistanceMetricType::New();
  typename DistanceMetricType::OriginType origin( ImageDimension );

  // Check for closely located seeds and delete
  unsigned int count = 0;
  it = allSeeds.begin();
  while ( it != allSeeds.end() )
    {
    flag = true;
    pt = it->second;
    for( unsigned int i = 0 ; i < ImageDimension; i++ )
      {
      mv[i] = pt[i];
      origin[i] = pt[i];
      }

    m_Tree->Search( mv, numberOfNeighbors, neighbors ) ;
    for ( unsigned int i = 0 ; i < neighbors.size() ; ++i )
      {
      distanceMetric->SetOrigin( origin );
      dist = distanceMetric->Evaluate( m_Tree->GetMeasurementVector( neighbors[i]) );
      if ( ( dist < m_LargestCellRadius ) && ( neighbors[i] < count ) )
        {
        flag = false;
        }
      }
    if ( flag )
      {
      seeds[it->first] = it->second;
      RemoveRegion( it->second );
      }

    ++it;
    count++;
    }
}


template < class TFeatureImage, class TInputImage, class TSegmentImage >
void
SeedExtraction< TFeatureImage, TInputImage, TSegmentImage >::
GenerateData()
{
  m_Spacing = this->GetInput()->GetSpacing();
  m_Size = this->GetInput()->GetLargestPossibleRegion().GetSize();

  // Run a connected component analysis on the foreground image
  // Eliminate components below a certain size threshold
  unsigned int NumOfLabels;
  {
    ConnectedComponentFilterPointer m_labelFilter =
      ConnectedComponentFilterType::New();
    m_labelFilter->SetInput ( m_fgMap );
    m_labelFilter->SetBackgroundValue ( 0 );
    m_labelFilter->Update();

    RelabelComponentFilterPointer m_relabelFilter =
      RelabelComponentFilterType::New();
    m_relabelFilter->SetInput ( m_labelFilter->GetOutput() );
    m_relabelFilter->SetMinimumObjectSize ( 100 );
    m_relabelFilter->Update();
    NumOfLabels = m_relabelFilter->GetNumberOfObjects();

    m_relabelfgMap = m_relabelFilter->GetOutput();
    m_relabelfgMap->DisconnectPipeline();
  }

  InputCastPointer caster = InputCastType::New();
  caster->SetInput( this->GetInput() );
  caster->Update();
  m_Input = caster->GetOutput();

  this->ImageSeeds();
  this->NonMaximalSeeds();

  // Identify the seeded components and store in a vector
  // Components[seedId] = 0 if it is seeded
  // The idea is to identify all possible unseeded components and add it later
  vnlVectorType Components;
  Components.set_size ( NumOfLabels + 1 );
  for ( unsigned int label = 1; label <= NumOfLabels; label++ )
    Components[label] = static_cast<float> ( label );
  Components[0] = 0;

  SeedIteratorType loc;
  ImagePointType pt;
  ImageIndexType idx;
  unsigned int val;

  loc = seeds.end();
  while ( loc != seeds.begin() )
  {
    --loc;
    pt = loc->second;
    m_relabelfgMap->TransformPhysicalPointToIndex( pt, idx );
    val = m_relabelfgMap->GetPixel ( idx );
    Components[val] = 0;
  }

  SegmentIteratorType fIt ( m_relabelfgMap, m_relabelfgMap->GetLargestPossibleRegion() );
  for ( fIt.GoToBegin(); !fIt.IsAtEnd(); ++fIt )
    if ( fIt.Get() )
    {
      fIt.Set ( Components[fIt.Get() ] );
    }

  // Run connected components to serialize the output components
  RelabelComponentFilterPointer m_relabelFilter =
    RelabelComponentFilterType::New();
  m_relabelFilter->SetInput ( m_relabelfgMap );
  m_relabelFilter->Update();

  this->GraftOutput ( m_relabelFilter->GetOutput() );
}


template < class TFeatureImage, class TInputImage, class TSegmentImage >
void
SeedExtraction< TFeatureImage, TInputImage, TSegmentImage >::
PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os,indent );
  os << indent << "Class Name:          " << GetNameOfClass( ) << std::endl;
  os << indent << "Largest Cell Radius: " << GetLargestCellRadius() <<
    std::endl;
}

} /* end namespace itk */

#endif

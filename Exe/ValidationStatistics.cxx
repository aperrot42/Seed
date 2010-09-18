#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"
#include "itkSmartPointer.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"

#include <iostream>
//#include <sstream>


int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << ": GroundTruthImage inputSeeds(.txt:position value) outputStats(.txt) pixelcoord-worldcoord(0-1)"
              << std::endl;
    return EXIT_FAILURE;
    }

    // Statistics Typedef
  typedef unsigned int NucleiIdType;
  typedef unsigned int NucleiCountType;
  typedef unsigned int SeedCountType;

  // Define the dimension of the images
  const unsigned int Dimension = 3;

  // Declare the types of the images
  //input
  typedef unsigned int      SegmentPixelType;
  typedef itk::Image< SegmentPixelType, Dimension>  SegmentImageType;
  
  // input image reader
  typedef itk::ImageFileReader< SegmentImageType > SegmentReaderType;
  
  // maxima filter
  typedef itk::MinimumMaximumImageCalculator<SegmentImageType> MaximaFilterType;

  typedef itk::ScalarConnectedComponentImageFilter< SegmentImageType, SegmentImageType >
    ConnectedComponentFilterType;

  //reading input image
  std::cout << "reading ground truth label image" << std::endl;
  SegmentReaderType::Pointer reader = SegmentReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();
  SegmentImageType::Pointer groundTruthImage = reader->GetOutput();


  ConnectedComponentFilterType::Pointer connectedComponentsFilter =
    ConnectedComponentFilterType::New();

  // make sure about internal storage type
  connectedComponentsFilter->SetInput(reader->GetOutput());
  connectedComponentsFilter->Update();
  connectedComponentsFilter->GetObjectCount();
  //displaying informations about the input image
  std::cout << "Spacing of segmented image" 
            << reader->GetOutput()->GetSpacing() << std::endl;
  std::cout << "Origin of segmented image"
            << reader->GetOutput()->GetOrigin() << std::endl;

  // looking for the maxima in the label image
  MaximaFilterType::Pointer maximaFilter = MaximaFilterType::New();
  maximaFilter->SetImage(groundTruthImage);
  maximaFilter->ComputeMaximum();
  NucleiCountType numberOfNuclei = maximaFilter->GetMaximum();
  std::cout << " Number of Nuclei in the image : " << numberOfNuclei << std::endl;
  std::cout << " or " << connectedComponentsFilter->GetObjectCount() << std::endl;

  std::vector<SeedCountType> seedsPerNuclei;
  seedsPerNuclei.resize(numberOfNuclei);
  NucleiIdType nucleiId;
  SeedCountType seedsOutOfNuclei = 0;
  SeedCountType overdectectedNuclei = 0;
  SeedCountType undetectedNucleis = numberOfNuclei;
  

  
  
  std::cout << "Read Seeds from " << argv[2] << std::endl;
  std::string line;
  //SeedPixelType ConfidencePix;
  SegmentImageType::PointType ptSeed;
  SegmentImageType::IndexType indexSeed;

  float fConfidenceSeed;
  std::stringstream sStream;
  std::ifstream seedFile (argv[2]);
  if (seedFile.is_open())
  {
    while (! seedFile.eof() )
    {
      getline (seedFile,line);
      if (line != "") // last line is empty...
      {
      // reading text file
      sStream.clear();
      sStream.str("");
      sStream.str(line);
      sStream >> ptSeed[0]
              >> ptSeed[1]
              >> ptSeed[2]
              >> fConfidenceSeed;
      // get pixel index and test value of Ground truth at this point
      groundTruthImage->TransformPhysicalPointToIndex( ptSeed, indexSeed );
      nucleiId = groundTruthImage->GetPixel(indexSeed);

      // Statistics

      // if there is a cell at this point
      if ( nucleiId )
         {
         
         // 0 is background in the image, whereas in the vector, it is the first cell
         if (seedsPerNuclei[nucleiId-1] == 0)
           {
           // we detected a new nuclei : good
           --undetectedNucleis;
           }
         else
           {
           if (seedsPerNuclei[nucleiId-1] == 1)
             {
             // if we already detected that one once,
             // we have an overdectectedNuclei
             ++overdectectedNuclei;
             }
           }
         seedsPerNuclei.at(nucleiId-1) = seedsPerNuclei[nucleiId-1] + 1;
         //++(seedsPerNuclei.at(nucleiId-1));
         }
       else
         {
         // we are in the choux
         ++seedsOutOfNuclei;
         }

      }
    }
    seedFile.close();
  }
  else
    {
    std::cerr << "Unable to open file";
    return EXIT_FAILURE;
    }

std::cout << "\t" << "| correct "
          << "\t" << "| outside "
          << "\t" << "| overdetect "
          << "\t" << "| undetected "
          << "\t" << "| total amount of nuclei |"
          <<std::endl;
std::cout << "\t" << "| "
          << "\t" << numberOfNuclei - (undetectedNucleis + overdectectedNuclei)
          << "\t" << "| "
          << "\t" << seedsOutOfNuclei
          << "\t" << "| "
          << "\t" << overdectectedNuclei
          << "\t" << "| "
          << "\t" << undetectedNucleis
          << "\t" << "| "
          << "\t" << numberOfNuclei
          << "\t" << "| "
          <<std::endl;

  std::ofstream statFile (argv[3]);
  if (statFile.is_open())
    {
    // we save the name of file for the seeds
    statFile << argv[2] << std::endl;
    
    statFile << "\t" << "| correct "
          << "\t" << "| outside "
          << "\t" << "| overdetect "
          << "\t" << "| undetected "
          << "\t" << "| total amount of nuclei |"
          <<std::endl;
    // we save the stats
    statFile << "\t" << "| "
          << "\t" << numberOfNuclei - (undetectedNucleis + overdectectedNuclei)
          << "\t" << "| "
          << "\t" << seedsOutOfNuclei
          << "\t" << "| "
          << "\t" << overdectectedNuclei
          << "\t" << "| "
          << "\t" << undetectedNucleis
          << "\t" << "| "
          << "\t" << numberOfNuclei
          << "\t" << "| "
          <<std::endl;
    statFile.close();
    }
  else
    {
    std::cerr << "Unable to write file";
    return EXIT_FAILURE;
    }
    
  return EXIT_SUCCESS;
}


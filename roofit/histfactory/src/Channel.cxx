// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, George Lewis
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/** \class RooStats::HistFactory::Channel
 *  \ingroup HistFactory
  This class encapsulates all information for the statistical interpretation of one experiment.
  It can be combined with other channels (e.g. for the combination of multiple experiments, or
  to constrain nuisance parameters with information obtained in a control region).
  A channel contains one or more samples which describe the contribution from different processes
  to this measurement.
*/



#include "RooStats/HistFactory/Channel.h"
#include "HFMsgService.h"
#include <cstdlib>

#include "TFile.h"
#include "TKey.h"
#include "TTimeStamp.h"

#include "RooStats/HistFactory/HistFactoryException.h"

using std::ofstream;

RooStats::HistFactory::Channel::Channel(std::string ChanName, std::string ChanInputFile) :
  fName( ChanName ), fInputFile( ChanInputFile )
{
  // create channel with given name and input file
}

namespace RooStats{
  namespace HistFactory{
    //BadChannel = Channel();
    Channel BadChannel;
    //    BadChannel.Name = "BadChannel"; // = Channel(); //.Name = "BadChannel";
  }
}


void RooStats::HistFactory::Channel::AddSample( RooStats::HistFactory::Sample sample )
{
  // add fully configured sample to channel

  sample.SetChannelName( GetName() );
  fSamples.push_back( sample );
}

void RooStats::HistFactory::Channel::Print( std::ostream& stream ) {
  // print information of channel to given stream

  stream << "\t Channel Name: " << fName
    << "\t InputFile: " << fInputFile
    << std::endl;

  stream << "\t Data:" << std::endl;
  fData.Print( stream );


  stream << "\t statErrorConfig:" << std::endl;
  fStatErrorConfig.Print( stream );


  if( !fSamples.empty() ) {

    stream << "\t Samples: " << std::endl;
    for( unsigned int i = 0; i < fSamples.size(); ++i ) {
      fSamples.at(i).Print( stream );
    }
  }


  stream << "\t End of Channel " << fName <<  std::endl;


}


void RooStats::HistFactory::Channel::PrintXML( std::string const& directory, std::string const& prefix ) const {

  // Create an XML file for this channel
  cxcoutPHF << "Printing XML Files for channel: " << GetName() << std::endl;

  std::string XMLName = prefix + fName + ".xml";
  if(!directory.empty()) XMLName = directory + "/" + XMLName;

  ofstream xml( XMLName.c_str() );

  // Add the time
  xml << "<!--" << std::endl;
  xml << "This xml file created automatically on: " << std::endl;
  // LM: use TTimeStamp since time_t does not work on Windows
  TTimeStamp t;
  UInt_t year = 0;
  UInt_t month = 0;
  UInt_t day = 0;
  t.GetDate(true, 0, &year, &month, &day);
  xml << year << '-'
      << month << '-'
      << day
      << std::endl;
  xml << "-->" << std::endl;

  // Add the DOCTYPE
  xml << "<!DOCTYPE Channel  SYSTEM 'HistFactorySchema.dtd'>  " << std::endl << std::endl;

  // Add the Channel
  xml << "  <Channel Name=\"" << fName << "\" InputFile=\"" << fInputFile << "\" >" << std::endl << std::endl;

  fData.PrintXML( xml );
  for(auto const& data : fAdditionalData) {
    data.PrintXML(xml);
  }

  fStatErrorConfig.PrintXML( xml );
  /*
  xml << "    <StatErrorConfig RelErrorThreshold=\"" << fStatErrorConfig.GetRelErrorThreshold() << "\" "
      << "ConstraintType=\"" << Constraint::Name( fStatErrorConfig.GetConstraintType() ) << "\" "
      << "/> " << std::endl << std::endl;
  */

  for(auto const& sample : fSamples) {
    sample.PrintXML( xml );
    xml << std::endl << std::endl;
  }

  xml << std::endl;
  xml << "  </Channel>  " << std::endl;
  xml.close();

  cxcoutPHF << "Finished printing XML files" << std::endl;

}


void RooStats::HistFactory::Channel::SetData( std::string DataHistoName, std::string DataInputFile, std::string DataHistoPath ) {
  // set data for this channel by specifying the name of the histogram,
  // the external ROOT file and the path to the histogram inside the ROOT file

  fData.SetHistoName( DataHistoName );
  fData.SetInputFile( DataInputFile );
  fData.SetHistoPath( DataHistoPath );

}



void RooStats::HistFactory::Channel::SetData( TH1* hData ) {
  // set data directly to some histogram
  fData.SetHisto( hData );
}

void RooStats::HistFactory::Channel::SetData( double val ) {

  // For a NumberCounting measurement only
  // Set the value of data in a particular channel
  //
  // Internally, this simply creates a 1-bin TH1F for you

  std::string DataHistName = fName + "_data";

  // Histogram has 1-bin (hard-coded)
  TH1F* hData = new TH1F( DataHistName.c_str(), DataHistName.c_str(), 1, 0, 1 );
  hData->SetBinContent( 1, val );

  // Set the histogram of the internally held data
  // node of this channel to this newly created histogram
  SetData( hData );

}


void RooStats::HistFactory::Channel::SetStatErrorConfig( double StatRelErrorThreshold, Constraint::Type StatConstraintType ) {

  fStatErrorConfig.SetRelErrorThreshold( StatRelErrorThreshold );
  fStatErrorConfig.SetConstraintType( StatConstraintType );

}

void RooStats::HistFactory::Channel::SetStatErrorConfig( double StatRelErrorThreshold, std::string StatConstraintType ) {

  fStatErrorConfig.SetRelErrorThreshold( StatRelErrorThreshold );
  fStatErrorConfig.SetConstraintType( Constraint::GetType(StatConstraintType) );

}

bool isInteger(float number)
{
    return std::abs(number - std::round(number)) < std::numeric_limits<float>::epsilon();
}


TH1* RooStats::HistFactory::Channel::Rebin(TH1* original)
{
    if (!fRebinning) return original;

    float left_edge = fHistoBinLow;
    float right_edge = fHistoBinHigh;

    if (std::isnan(left_edge))
    {
      left_edge = original->GetXaxis()->GetXmin();
    }

    if (std::isnan(right_edge))
    {
      right_edge = original->GetXaxis()->GetXmax();
    }
    
    float original_bin_width = (original->GetXaxis()->GetBinWidth(1));
    int assert_check_left =  ((left_edge - original->GetXaxis()->GetXmin())) / original_bin_width;
    int assert_check_right = (original->GetXaxis()->GetXmax() - right_edge) / original_bin_width;

    std::cout << "fHistBinLow: " << left_edge << std::endl;
    std::cout << "fHistBinHigh: " << right_edge << std::endl;
    std::cout << "original_bin_width: " << original_bin_width << std::endl;

    if (!isInteger(assert_check_left) || !isInteger(assert_check_right))
    {
        std::cout << "Error: The left_edge and right_edge are not multiples of the original bin width" << std::endl;
        return original;
    }

    int number_of_remove_bins = ((left_edge - original->GetXaxis()->GetXmin()) + (original->GetXaxis()->GetXmax() - right_edge)) / original_bin_width;
    std::cout << "number of remove bins: " << number_of_remove_bins << std::endl;
    int new_nbins = (original->GetNbinsX() - number_of_remove_bins) / fRebin;
    std::cout << "new_nbins: " << new_nbins << std::endl;
    TH1F* original_new = new TH1F(original->GetName(), original->GetTitle(), new_nbins, left_edge, right_edge);
    int skipped_bins_left = (left_edge - original->GetXaxis()->GetXmin()) / original_bin_width;
    for (int i = 1; i <= new_nbins; ++i)
    {
        int bin = i + skipped_bins_left;
        original_new->SetBinContent(i, original->GetBinContent(bin));
        original_new->SetBinError(i, original->GetBinError(bin));
    }

    return original_new;
}



void RooStats::HistFactory::Channel::CollectHistograms() {

  // Loop through all Samples and Systematics
  // and collect all necessary histograms

  // Handles to open files for collecting histograms
  std::map<std::string,std::unique_ptr<TFile>> fileHandles;

  // Get the Data Histogram:

  if( !fData.GetInputFile().empty() ) {
    fData.SetHisto( Rebin(GetHistogram(fData.GetInputFile(),
             fData.GetHistoPath(),
             fData.GetHistoName(),
             fileHandles) ));
  }


  // Collect any histograms for additional Datasets
  for(auto& data : fAdditionalData) {
    if( !data.GetInputFile().empty() ) {
      data.SetHisto( Rebin(GetHistogram(data.GetInputFile(), data.GetHistoPath(), data.GetHistoName(), fileHandles) ));
    }
  }

  // Get the histograms for the samples:
  for( unsigned int sampItr = 0; sampItr < fSamples.size(); ++sampItr ) {

    RooStats::HistFactory::Sample& sample = fSamples.at( sampItr );


    // Get the nominal histogram:
    cxcoutDHF << "Collecting Nominal Histogram" << std::endl;
    TH1* Nominal =  Rebin(GetHistogram(sample.GetInputFile(),
             sample.GetHistoPath(),
             sample.GetHistoName(),
             fileHandles));

    sample.SetHisto( Nominal );

    // printf("number of histoSys: %d\n", sample.GetHistoSysList().size());
    
    for( unsigned int histoSysItr = 0; histoSysItr < sample.GetHistoSysList().size(); ++histoSysItr ) {
          RooStats::HistFactory::HistoSys& histoSys = sample.GetHistoSysList().at( histoSysItr );
            if (histoSys.GetSymmetrize())
            {

              
              // TODO: implement calculation of symmetrized histogram
              if (!histoSys.GetNormPlusShape())
              {
                TFile* file = TFile::Open(histoSys.GetInputFileHigh().c_str(), "UPDATE");
                TH1* h_top = Rebin(GetHistogram(histoSys.GetInputFileHigh(),
                  histoSys.GetHistoPathHigh(),
                  histoSys.GetHistoNameHigh(),
                  fileHandles));

                TH1* histNominal = Rebin(GetHistogram(sample.GetInputFile(),
                  sample.GetHistoPath(),
                  sample.GetHistoName(),
                  fileHandles));

                auto* h_new = (TH1*)h_top->Clone();
                h_new->Add(histNominal, -1);
                TH1* h_temp = (TH1*)histNominal->Clone(histoSys.GetHistoNameLow().c_str());
                // h_temp->SetName(histoSys.GetHistoNameLow().c_str());
                h_temp->Add(h_new, -1);
                file->cd(histoSys.GetHistoPathLow().c_str());
                bool write_histogram = true;
                if (write_histogram)  h_temp->Write();
                // printf("Setting histogram\n");
                histoSys.SetHistoLow(h_temp);
                // printf("Set histogram\n");

                file->Close();
              }
              else
              {
                TFile* file = TFile::Open(histoSys.GetInputFileHigh().c_str(), "UPDATE");
                TH1* h_top = Rebin(GetHistogram(histoSys.GetInputFileHigh(),
                  histoSys.GetHistoPathHigh(),
                  histoSys.GetHistoNameHigh(),
                  fileHandles));

                TH1* histNominal = Rebin(GetHistogram(sample.GetInputFile(),
                  sample.GetHistoPath(),
                  sample.GetHistoName(),
                  fileHandles));

                double norm_factor = h_top->Integral()/histNominal->Integral();

                auto* h_new = (TH1*)h_top->Clone(histoSys.GetHistoNameHigh().c_str());
                h_new->Scale(1 / norm_factor);

                // h_new->Add(histNominal, -1);
                TH1* h_temp = (TH1*)histNominal->Clone(histoSys.GetHistoNameLow().c_str());
                // h_temp->SetName(histoSys.GetHistoNameLow().c_str());
                h_temp->Scale(2);
                h_temp->Add(h_new, -1);
                // h_new->Scale(norm_factor*norm_factor);
                // double lo_scale = 2 - norm_factor;
                // h_temp->Scale(lo_scale);
                file->cd(histoSys.GetHistoPathLow().c_str());

                // trying to remove histogram
                auto* h1 = (TH1*)file->Get(histoSys.GetHistoNameLow().c_str());
                if (h1) h1->SetDirectory(nullptr);
                auto* h2 = (TH1*)file->Get(histoSys.GetHistoNameHigh().c_str());
                if (h2) h2->SetDirectory(nullptr);

                bool write_histogram = true;
                if (write_histogram)  
                  {
                    h_temp->Write();
                    h_new->Write();
                  }
                // printf("Setting histogram\n");
                histoSys.SetHistoLow(h_temp);
                histoSys.SetHistoHigh(h_new);
                // printf("Set histogram\n");

                // histoSys.Set

                RooStats::HistFactory::OverallSys sys;
                sys.SetName(histoSys.GetName());
                double sys_low = 2 - norm_factor;
                double sys_high = norm_factor;
                sys.SetLow(sys_low);
                sys.SetHigh(sys_high);
                sample.GetOverallSysList().push_back(sys);


                file->Close();
              }
            }
    }

    // Get the StatError Histogram (if necessary)
    if( sample.GetStatError().GetUseHisto() ) {
      sample.GetStatError().SetErrorHist( Rebin(GetHistogram(sample.GetStatError().GetInputFile(),
                         sample.GetStatError().GetHistoPath(),
                         sample.GetStatError().GetHistoName(),
                         fileHandles) ));
    }

    printf("number of histoSys: %d\n", sample.GetHistoSysList().size());

    // Get the HistoSys Variations:
    for( unsigned int histoSysItr = 0; histoSysItr < sample.GetHistoSysList().size(); ++histoSysItr ) 
    {
      RooStats::HistFactory::HistoSys& histoSys = sample.GetHistoSysList().at( histoSysItr );

      if (!histoSys.GetSymmetrize())
      {
        histoSys.SetHistoLow( Rebin(GetHistogram(histoSys.GetInputFileLow(),
                  histoSys.GetHistoPathLow(),
                  histoSys.GetHistoNameLow(),
                  fileHandles) ));
      // }

        histoSys.SetHistoHigh( Rebin(GetHistogram(histoSys.GetInputFileHigh(),
                  histoSys.GetHistoPathHigh(),
                  histoSys.GetHistoNameHigh(),
                  fileHandles) ));
      }

      if (histoSys.GetNormPlusShape() and (!histoSys.GetSymmetrize()))
      {
        RooStats::HistFactory::OverallSys sys;
        sys.SetName(histoSys.GetName());
        double sys_low = histoSys.GetHistoLow()->Integral()/sample.GetHisto()->Integral();
        double sys_high = histoSys.GetHistoHigh()->Integral()/sample.GetHisto()->Integral();
        // double sys_high = 2 - sys_low;
        sys.SetLow(sys_low);
        sys.SetHigh(sys_high);
        sample.GetOverallSysList().push_back(sys);

        TFile* file = TFile::Open(histoSys.GetInputFileHigh().c_str(), "UPDATE");

        file->cd(histoSys.GetHistoPathLow().c_str()); // TODO: do this for both of files

        auto* histo_new_high = (TH1*)histoSys.GetHistoHigh()->Clone(histoSys.GetHistoNameHigh().c_str());
        histo_new_high->Scale(1 / sys_high);
        histoSys.SetHistoHigh(histo_new_high);


        auto* histo_new_low = (TH1*)histoSys.GetHistoLow()->Clone(histoSys.GetHistoNameLow().c_str());
        histo_new_low->Scale(1 / sys_low);
        histoSys.SetHistoLow(histo_new_low);

        bool write_histogram = true;
        if (write_histogram)  
        {
          histo_new_high->Write();
          histo_new_low->Write();
        }


        file->Close();

        // std::string shape_name = "Test_name";
        // RooStats::HistFactory::ShapeFactor factor;
        // factor.SetName(shape_name);
        // sample.GetShapeFactorList().push_back(factor);
      }
    } // End Loop over HistoSys


      // Get the HistoFactor Variations:
    for( unsigned int histoFactorItr = 0; histoFactorItr < sample.GetHistoFactorList().size(); ++histoFactorItr ) {

      RooStats::HistFactory::HistoFactor& histoFactor = sample.GetHistoFactorList().at( histoFactorItr );

      histoFactor.SetHistoLow( Rebin(GetHistogram(histoFactor.GetInputFileLow(),
                   histoFactor.GetHistoPathLow(),
                   histoFactor.GetHistoNameLow(),
                   fileHandles) ));

      histoFactor.SetHistoHigh( Rebin(GetHistogram(histoFactor.GetInputFileHigh(),
                    histoFactor.GetHistoPathHigh(),
                    histoFactor.GetHistoNameHigh(),
                    fileHandles) ));
    } // End Loop over HistoFactor


      // Get the ShapeSys Variations:
    for( unsigned int shapeSysItr = 0; shapeSysItr < sample.GetShapeSysList().size(); ++shapeSysItr ) {

      RooStats::HistFactory::ShapeSys& shapeSys = sample.GetShapeSysList().at( shapeSysItr );

      shapeSys.SetErrorHist( Rebin(GetHistogram(shapeSys.GetInputFile(),
                 shapeSys.GetHistoPath(),
                 shapeSys.GetHistoName(),
                 fileHandles) ));
    } // End Loop over ShapeSys


    // Get any initial shape for a ShapeFactor
    for( unsigned int shapeFactorItr = 0; shapeFactorItr < sample.GetShapeFactorList().size(); ++shapeFactorItr ) {

      RooStats::HistFactory::ShapeFactor& shapeFactor = sample.GetShapeFactorList().at( shapeFactorItr );

      // Check if we need an InitialShape
      if( shapeFactor.HasInitialShape() ) {
   TH1* hist = Rebin(GetHistogram( shapeFactor.GetInputFile(), shapeFactor.GetHistoPath(),
              shapeFactor.GetHistoName(), fileHandles ));
   shapeFactor.SetInitialShape( hist );
      }

    } // End Loop over ShapeFactor

  } // End Loop over Samples
}


bool RooStats::HistFactory::Channel::CheckHistograms() const {

  // Check that all internal histogram pointers
  // are properly configured (ie that they're not nullptr)

    if( fData.GetHisto() == nullptr && !fData.GetInputFile().empty() ) {
      cxcoutEHF << "Error: Data Histogram for channel " << GetName() << " is nullptr." << std::endl;
      return false;
    }

    // Get the histograms for the samples:
    for( unsigned int sampItr = 0; sampItr < fSamples.size(); ++sampItr ) {

      const RooStats::HistFactory::Sample& sample = fSamples.at( sampItr );

      // Get the nominal histogram:
      if( sample.GetHisto() == nullptr ) {
   cxcoutEHF << "Error: Nominal Histogram for sample " << sample.GetName() << " is nullptr." << std::endl;
   return false;
      }
      else {

   // Check if any bins are negative
   std::vector<int> NegativeBinNumber;
   std::vector<double> NegativeBinContent;
   const TH1* histNominal = sample.GetHisto();
   for(int ibin=1; ibin<=histNominal->GetNbinsX(); ++ibin) {
     if(histNominal->GetBinContent(ibin) < 0) {
       NegativeBinNumber.push_back(ibin);
       NegativeBinContent.push_back(histNominal->GetBinContent(ibin));
     }
   }
   if(!NegativeBinNumber.empty()) {
     cxcoutWHF << "WARNING: Nominal Histogram " << histNominal->GetName() << " for Sample = " << sample.GetName()
          << " in Channel = " << GetName() << " has negative entries in bin numbers = ";

     for(unsigned int ibin=0; ibin<NegativeBinNumber.size(); ++ibin) {
       if(ibin>0) std::cout << " , " ;
       std::cout << NegativeBinNumber[ibin] << " : " << NegativeBinContent[ibin] ;
     }
     std::cout << std::endl;
   }

      }

      // Get the StatError Histogram (if necessary)
      if( sample.GetStatError().GetUseHisto() ) {
   if( sample.GetStatError().GetErrorHist() == nullptr ) {
     cxcoutEHF << "Error: Statistical Error Histogram for sample " << sample.GetName() << " is nullptr." << std::endl;
     return false;
   }
      }


      // Get the HistoSys Variations:
      for( unsigned int histoSysItr = 0; histoSysItr < sample.GetHistoSysList().size(); ++histoSysItr ) {

        const RooStats::HistFactory::HistoSys& histoSys = sample.GetHistoSysList().at( histoSysItr );

          // if (histoSys.GetSymmetrize())
          // {
          //   // TODO: implement calculation of symmetrized histogram

          //   TFile* file = TFile::Open(histoSys.GetInputFileHigh().c_str(), "UPDATE");
            
          //   TH1* h_top = dynamic_cast<TH1*>(file->Get((histoSys.GetHistoPathHigh()+ "/" + histoSys.GetHistoNameHigh()).c_str()));
          //   const TH1* histNominal = sample.GetHisto();
          //   auto* h_new = (TH1*)h_top->Clone();
          //   h_new->Add(histNominal, -1);
          //   auto* h_temp = (TH1*)histNominal->Clone(histoSys.GetHistoNameLow().c_str());
          //   h_temp->Add(h_new, -1);

          //   file->cd(histoSys.GetHistoPathLow().c_str());
          //   h_temp->Write();

          //   file->Close();
          // }

        if( histoSys.GetHistoLow() == nullptr ) {
          cxcoutEHF << "Error: HistoSyst Low for Systematic " << histoSys.GetName()
                << " in sample " << sample.GetName() << " is nullptr." << std::endl;
          return false;
        }
        if( histoSys.GetHistoHigh() == nullptr ) {
          cxcoutEHF << "Error: HistoSyst High for Systematic " << histoSys.GetName()
                << " in sample " << sample.GetName() << " is nullptr." << std::endl;
          return false;
        }

      } // End Loop over HistoSys


      // Get the HistoFactor Variations:
      for( unsigned int histoFactorItr = 0; histoFactorItr < sample.GetHistoFactorList().size(); ++histoFactorItr ) {

   const RooStats::HistFactory::HistoFactor& histoFactor = sample.GetHistoFactorList().at( histoFactorItr );

   if( histoFactor.GetHistoLow() == nullptr ) {
     cxcoutEHF << "Error: HistoSyst Low for Systematic " << histoFactor.GetName()
          << " in sample " << sample.GetName() << " is nullptr." << std::endl;
     return false;
   }
   if( histoFactor.GetHistoHigh() == nullptr ) {
     cxcoutEHF << "Error: HistoSyst High for Systematic " << histoFactor.GetName()
          << " in sample " << sample.GetName() << " is nullptr." << std::endl;
     return false;
   }

      } // End Loop over HistoFactor


      // Get the ShapeSys Variations:
      for( unsigned int shapeSysItr = 0; shapeSysItr < sample.GetShapeSysList().size(); ++shapeSysItr ) {

   const RooStats::HistFactory::ShapeSys& shapeSys = sample.GetShapeSysList().at( shapeSysItr );

   if( shapeSys.GetErrorHist() == nullptr ) {
     cxcoutEHF << "Error: HistoSyst High for Systematic " << shapeSys.GetName()
          << " in sample " << sample.GetName() << " is nullptr." << std::endl;
     return false;
   }

      } // End Loop over ShapeSys

    } // End Loop over Samples

  return true;
}



/// Open a file and copy a histogram
/// \param InputFile File where the histogram resides.
/// \param HistoPath Path of the histogram in the file.
/// \param HistoName Name of the histogram to retrieve.
/// \param lsof List of open files. Helps to prevent opening and closing a file hundreds of times.
TH1* RooStats::HistFactory::Channel::GetHistogram(std::string InputFile, std::string HistoPath, std::string HistoName, std::map<std::string,std::unique_ptr<TFile>>& lsof) {

  cxcoutPHF << "Getting histogram " << InputFile << ":" << HistoPath << "/" << HistoName << std::endl;

  auto& inFile = lsof[InputFile];
  if (!inFile || !inFile->IsOpen()) {
    inFile.reset( TFile::Open(InputFile.c_str()) );
    if ( !inFile || !inFile->IsOpen() ) {
      cxcoutEHF << "Error: Unable to open input file: " << InputFile << std::endl;
      throw hf_exc();
    }
    cxcoutIHF << "Opened input file: " << InputFile << ": " << std::endl;
  }

  TDirectory* dir = inFile->GetDirectory(HistoPath.c_str());
  if (dir == nullptr) {
    cxcoutEHF << "Histogram path '" << HistoPath
        << "' wasn't found in file '" << InputFile << "'." << std::endl;
    throw hf_exc();
  }

  // Have to read histograms via keys, to ensure that the latest-greatest
  // name cycle is read from file. Otherwise, they might come from memory.
  auto key = dir->GetKey(HistoName.c_str());
  if (key == nullptr) {
    cxcoutEHF << "Histogram '" << HistoName
        << "' wasn't found in file '" << InputFile
        << "' in directory '" << HistoPath << "'." << std::endl;
    throw hf_exc();
  }

  auto hist = key->ReadObject<TH1>();
  if( !hist ) {
    cxcoutEHF << "Histogram '" << HistoName
        << "' wasn't found in file '" << InputFile
        << "' in directory '" << HistoPath << "'." << std::endl;
    throw hf_exc();
  }


  TH1 * ptr = static_cast<TH1 *>(hist->Clone());

  if(!ptr){
    std::cerr << "Not all necessary info are set to access the input file. Check your config" << std::endl;
    std::cerr << "filename: " << InputFile
         << "path: " << HistoPath
         << "obj: " << HistoName << std::endl;
    throw hf_exc();
  }

  ptr->SetDirectory(nullptr);


#ifdef DEBUG
  std::cout << "Found Histogram: " << HistoName " at address: " << ptr
       << " with integral "   << ptr->Integral() << " and mean " << ptr->GetMean()
       << std::endl;
#endif

  // Done
  return ptr;

}

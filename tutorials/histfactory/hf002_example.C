/// \file
/// \ingroup tutorial_histfactory
/// A ROOT script demonstrating  an example of writing a HistFactory  model using c++ only.
///
/// \macro_code
/// \macro_output
///
/// \author George Lewis


#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/ModelConfig.h"
#include "TFile.h"
#include "TROOT.h"

using namespace RooStats;
using namespace HistFactory;

void hf002_example() {


  std::string InputFile = "./hf002_input_data.root";
  // in case the file is not found
  bool bfile = gSystem->AccessPathName(InputFile.c_str());
  if (bfile) {
     std::cout << "Input file is not found - run prepareHistFactory script " << std::endl;
     gROOT->ProcessLine(".! prepareHistFactory .");
     bfile = gSystem->AccessPathName(InputFile.c_str());
     if (bfile) {
        std::cout << "Still no " << InputFile << ", giving up.\n";
        exit(1);
     }
  }

  // Create the measurement
  Measurement meas("meas", "meas");

  meas.SetOutputFilePrefix( "./results/example_UsingC" );
  meas.SetPOI( "SigXsecOverSM" );
  meas.AddConstantParam("Lumi");

  meas.SetLumi( 1.0 );
  meas.SetLumiRelErr( 0.0 );
  meas.SetExportOnly( false );
  meas.SetBinHigh( 2 );

  // Create a channel

  Channel chan( "channel1" );
  chan.SetData( "data", InputFile );
  chan.SetStatErrorConfig( 0.05, "Poisson" );
  chan.SetGlobalOverallSys("syst1", 0.95, 1.05);


  // Now, create some samples


  // Create the signal sample
  Sample signal( "signal", "signal", InputFile );
  signal.ActivateStatError();
  signal.AddNormFactor("SigXsecOverSM", 1, 0, 3 );
  chan.AddSample( signal );

  // Background 1
  Sample background1("background1", "background1", InputFile );
  // Input file name can be skipped, in such case same input file as for nominal histogram will be used
  background1.AddHistoSys("weighted_modeling",
                            "weighted_modeling_up", "", "",
                              "weighted_modeling_down", "", "");
  background1.AddHistoSys("modeling",
                            "modeling_up", "", "", 
                              "modeling_down", "", "");
  background1.ActivateStatError();
  chan.AddSample( background1 );



  // Done with this channel
  // Add it to the measurement:
  meas.AddChannel( chan );

  // Collect the histograms from their files,
  meas.CollectHistograms();

  // Now, do the measurement and create RooWorkspace
  std::unique_ptr<RooWorkspace> ws{MakeModelAndMeasurementFast(meas)};

  RooStats::ModelConfig *modelConfig = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

  RooAbsPdf *pdf = modelConfig->GetPdf();
  RooArgSet globalObservables{*modelConfig->GetGlobalObservables()};

  // Perform the fit
  using namespace RooFit;
  std::unique_ptr<RooFitResult> result{
      pdf->fitTo(*ws->data("obsData"), Save(), PrintLevel(-1), GlobalObservables(globalObservables))
  };

  result->Print();

  bool saveFitResults = true;

  if (saveFitResults)
  {
    std::unique_ptr<TFile> file = std::make_unique<TFile>("fitResults.root", "RECREATE");
    result->Write("fitResult");
    file->Close();
  }

}


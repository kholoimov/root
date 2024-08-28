///
/// code which show how histogram model, built by HistFactory can be drawn
/// extra features: prefit and postfit uncertainties calculation for each bin
///
///

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/ModelConfig.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TROOT.h"
#include "TH1.h"
#include "TColor.h"
#include "TMath.h"

using namespace RooStats;
using namespace HistFactory;

const std::vector<Int_t> colors 
{
    TColor::GetColor("#3F90DA"), 
    TColor::GetColor("#FFA90E"), 
    TColor::GetColor("#BD1F01"),
    TColor::GetColor("#94A4A2"),  
    TColor::GetColor("#832DB6"),  
    TColor::GetColor("#A96B59"),
    TColor::GetColor("#E76300"),  
    TColor::GetColor("#B9AC70"),  
    TColor::GetColor("#717581"),  
    TColor::GetColor("#92DADD"), 
};


void DrawHistogramUsingMeasurement(const Measurement& meas)
{
    const auto channels = meas.GetChannels();


    for (auto& channel : channels)
    {
        std::unique_ptr<TCanvas> cv = std::make_unique<TCanvas>(channel.GetName().c_str(), channel.GetName().c_str(), 800, 600);

        const TH1* data_histogram = channel.GetData().GetHisto();

        std::unique_ptr<THStack> sample_stack = std::make_unique<THStack>();
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.6, 0.8, 0.9, 0.9);

        int color_index = 0;
        for (auto& sample : channel.GetSamples())
        {
            TH1* sample_histogram = static_cast<TH1*>(sample.GetHisto()->Clone());
            
            sample_histogram->SetFillColor(colors[color_index]);
            color_index++;

            legend->AddEntry(sample_histogram, sample.GetName().c_str(), "f");

            sample_stack->Add(sample_histogram); // error here
        }

        sample_stack->Draw("hist");

        int maximum_bin_data = data_histogram->GetMaximum();
        int maximum_bin_stack = sample_stack->GetMaximum();

        sample_stack->SetMaximum(1.2 * ((maximum_bin_data > maximum_bin_stack) ? maximum_bin_data : maximum_bin_stack)); // free some space for legend


        const_cast<TH1*>(data_histogram)->SetMarkerStyle(8);
        const_cast<TH1*>(data_histogram)->Draw("same P");

        legend->AddEntry(data_histogram, "Pseudodata", "p");
        legend->Draw();
        
        cv->Draw();
        cv->SaveAs((std::string(channel.GetName()) + "_hist_model.png").c_str());
    }
}



void DrawHistogramUsingRooWorkspace(RooWorkspace* ws, Measurement& meas, RooFitResult* result)
{

    for (auto& channel : meas.GetChannels())
    {
        const std::string channel_name = std::string(channel.GetName());

        const TH1* data_histogram = channel.GetData().GetHisto();
        const int numBins = data_histogram->GetNbinsX();

        const int number_of_samples = channel.GetSamples().size();
        std::vector<std::string> sample_names;
        sample_names.reserve(number_of_samples); // we don't want to use dynamic allocation :)

        std::unique_ptr<TCanvas> cv = std::make_unique<TCanvas>(channel.GetName().c_str(), channel.GetName().c_str(), 800, 600);
        cv->Divide(2, 1);

        cv->cd(1);


        std::unique_ptr<THStack> sample_stack = std::make_unique<THStack>();
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.6, 0.8, 0.9, 0.9);

        int color_index = 0;
        
        std::vector<float> totalBinHeight(numBins, 0);

        for (auto& sample : channel.GetSamples())
        {
            TH1* sample_histogram = static_cast<TH1*>(sample.GetHisto()->Clone());
            
            sample_histogram->SetFillColor(colors[color_index]);
            color_index++;

            legend->AddEntry(sample_histogram, sample.GetName().c_str(), "f");

            sample_names.push_back(sample.GetName().c_str());

            for (int i_bin = 0; i_bin < numBins; i_bin++)
            {
                totalBinHeight[i_bin] += sample_histogram->GetBinContent(i_bin + 1);
            }


            sample_stack->Add(sample_histogram); // error here
        }

        sample_stack->Draw("hist");

        int maximum_bin_data = data_histogram->GetMaximum();
        int maximum_bin_stack = sample_stack->GetMaximum();

        float y_max_limit = 1.2 * ((maximum_bin_data > maximum_bin_stack) ? maximum_bin_data : maximum_bin_stack); 

        sample_stack->SetMaximum(y_max_limit); // free some space for legend


        const_cast<TH1*>(data_histogram)->SetMarkerStyle(8);
        const_cast<TH1*>(data_histogram)->Draw("same P");

        std::vector<float> yields_uncertainties(numBins); /// uncertainty in each bin

        auto channel_pdf = static_cast<RooRealSumPdf*>(ws->pdf((channel_name + "_model")));
        auto variable = ws->var("obs_x_" + channel_name);
        auto observables = RooArgSet(*variable);

        auto all_params = channel_pdf->getParameters(observables);
        const int number_of_params = all_params->size();

        const auto prefitParamValues = result->floatParsInit();


        for (int i_bin = 0; i_bin < numBins; i_bin++)
        {
            variable->setBin(i_bin);
            const float bin_width = variable->getBinWidth(i_bin);


            std::vector<float> uncertainties_summed(number_of_params);

            for (int i_sample = 0; i_sample < number_of_samples; i_sample++)
            {
                auto sample_yield = RooProduct("tmp", "tmp", {channel_pdf->funcList()[i_sample], channel_pdf->coefList()[i_sample]});

                int param_index = 0;
                for (auto par : *all_params)
                {
                    RooRealVar* real_par = dynamic_cast<RooRealVar*>(par);
                    RooRealVar* original_param = dynamic_cast<RooRealVar*>(prefitParamValues.find(par->GetName()));

                    if (original_param == nullptr)
                    {
                        continue;
                    }

                    const float original_param_value = real_par->getVal();
                    
                    const float cen_val = original_param->getVal();
                    const float par_err = original_param->getError();

                    real_par->setVal(cen_val + par_err);
                    float upper_var = sample_yield.getVal(observables); // should observables be inside function call

                    real_par->setVal(cen_val - par_err);
                    float down_var = sample_yield.getVal(observables);

                    real_par->setVal(original_param_value);
                    sample_yield.getVal();

                    uncertainties_summed[param_index] += (upper_var - down_var) / 2;
                    param_index++;
                }
            }

            /// for prefit total uncertainty can be calculated just as sum of squares

            float total_sum = 0;
            for (auto unc_val : uncertainties_summed)
            {
                total_sum += unc_val * unc_val;
            }

            yields_uncertainties[i_bin] = TMath::Sqrt(total_sum) * bin_width;
        }

        std::vector<TBox> error_boxes(numBins);

        float left_canvas_edge = data_histogram->GetBinLowEdge(1);
        float right_canvas_edge = data_histogram->GetBinLowEdge(numBins) + data_histogram->GetBinWidth(numBins);

        for (int i_bin = 1; i_bin < numBins + 1; i_bin++) /// ROOT Histogram have additional underflow and overflow bins, which should be handled manually
        {
            float leftEdge = data_histogram->GetBinLowEdge(i_bin);
            float binWidth = data_histogram->GetBinWidth(i_bin);
            float rightEdge = leftEdge + binWidth;

            float central_value = totalBinHeight[i_bin - 1];

            float lowerEdge = central_value - yields_uncertainties[i_bin - 1];
            float upperEdge = central_value + yields_uncertainties[i_bin - 1];

            if (upperEdge > y_max_limit) upperEdge = y_max_limit;
            
            error_boxes[i_bin - 1] = TBox(leftEdge, lowerEdge, rightEdge, upperEdge);
            error_boxes[i_bin - 1].SetFillStyle(3004);
            error_boxes[i_bin - 1].SetFillColor(kGray + 3);
            error_boxes[i_bin - 1].Draw("same");
        }

        legend->AddEntry(data_histogram, "Pseudodata", "p");
        legend->Draw();

        cv->cd(2);


        /// Getting postfit sample values and postfit uncertainties
        std::vector<float> yields(numBins, 0); /// sum of samples in each bin
        std::vector<std::vector<float>> sample_values(number_of_samples, std::vector<float>(numBins)); /// value for each sample for each bin
        std::vector<float> postfit_yields_uncertainties(numBins); /// uncertainty in each bin
        

        RooArgList paramList;
        for (auto par : *all_params)
        {
            auto rrvInAbsReal = static_cast<RooRealVar const*>(par);
            if (rrvInAbsReal->getError() <= std::abs(rrvInAbsReal->getVal()) * std::numeric_limits<double>::epsilon()) continue;
            paramList.add(*rrvInAbsReal);
        }

        TMatrixDSym covMatrix = result->reducedCovarianceMatrix(paramList);

        int total_events = channel_pdf->expectedEvents(&observables);
        for (int i_bin = 0; i_bin < numBins; i_bin++)
        {
            variable->setBin(i_bin);
            const float bin_width = variable->getBinWidth(i_bin);

            TVectorD uncertainties_summed(paramList.size());

            for (int i_sample = 0; i_sample < number_of_samples; i_sample++)
            {
                auto sample_yield = RooProduct("tmp", "tmp", {channel_pdf->funcList()[i_sample], channel_pdf->coefList()[i_sample]});

                yields[i_bin] += sample_yield.getVal() * bin_width;
                sample_values[i_sample][i_bin] = sample_yield.getVal() * bin_width;

                int param_index = 0;
                for (auto par : paramList)
                {
                    RooRealVar* real_par = dynamic_cast<RooRealVar*>(par);
                    
                    const float cen_val = real_par->getVal();
                    const float par_err = real_par->getError();

                    real_par->setVal(cen_val + par_err);
                    float upper_var = sample_yield.getVal(observables); // should observables be inside function call

                    real_par->setVal(cen_val - par_err);
                    float down_var = sample_yield.getVal(observables);

                    real_par->setVal(cen_val);
                    sample_yield.getVal();

                    uncertainties_summed[param_index] += (upper_var - down_var) / 2;
                    param_index++;
                }
            }

            /// for prefit total uncertainty can be calculated just as sum of squares

            float total_sum = 0;
            printf("aborting here");

            total_sum = uncertainties_summed * (covMatrix * uncertainties_summed);

            postfit_yields_uncertainties[i_bin] = TMath::Sqrt(total_sum) * bin_width;

            // postfit_yields_uncertainties[i_bin] = channel_pdf->getPropagatedError(*result, observables) * total_events * bin_width;
        }

        std::unique_ptr<THStack> postfit_sample_stack = std::make_unique<THStack>();
        std::unique_ptr<TLegend> postfit_legend = std::make_unique<TLegend>(0.6, 0.8, 0.9, 0.9);
        std::vector<TH1F> sample_histograms_vec(number_of_samples);

        color_index = 0;

        for (auto value_vec : sample_values)
        {
            for (auto int_value : value_vec)
            {
                printf("%f\t", int_value);
            }
            printf("\n");
        }
        
        for (int i_sample = 0; i_sample < number_of_samples; i_sample++)
        {
            sample_histograms_vec[i_sample] = TH1F((sample_names[i_sample] + "postfit").c_str(), (sample_names[i_sample] + "postfit").c_str(), numBins, left_canvas_edge, right_canvas_edge);

            for (int i_bin = 1; i_bin <= numBins; i_bin++)
            {
                sample_histograms_vec[i_sample].SetBinContent(i_bin, sample_values[i_sample][i_bin - 1]);
            }
            
            sample_histograms_vec[i_sample].SetFillColor(colors[color_index]);
            color_index++;

            postfit_legend->AddEntry(&sample_histograms_vec[i_sample], sample_names[i_sample].c_str(), "f");

            postfit_sample_stack->Add(&sample_histograms_vec[i_sample]); // error here
        }

        postfit_legend->AddEntry(data_histogram, "Pseudodata", "p");

        postfit_sample_stack->Draw("hist");
        const_cast<TH1*>(data_histogram)->Draw("same p");
        postfit_legend->Draw();

        maximum_bin_stack = postfit_sample_stack->GetMaximum();

        y_max_limit = 1.2 * ((maximum_bin_data > maximum_bin_stack) ? maximum_bin_data : maximum_bin_stack); 
        postfit_sample_stack->SetMaximum(y_max_limit); 

        std::vector<TBox> postfit_error_boxes(numBins);

        for (int i_bin = 1; i_bin < numBins + 1; i_bin++) /// ROOT Histogram have additional underflow and overflow bins, which should be handled manually
        {
            float leftEdge = data_histogram->GetBinLowEdge(i_bin);
            float binWidth = data_histogram->GetBinWidth(i_bin);
            float rightEdge = leftEdge + binWidth;

            float central_value = yields[i_bin - 1];

            float lowerEdge = central_value - postfit_yields_uncertainties[i_bin - 1];
            float upperEdge = central_value + postfit_yields_uncertainties[i_bin - 1];

            if (upperEdge > y_max_limit) upperEdge = y_max_limit;
            
            postfit_error_boxes[i_bin - 1] = TBox(leftEdge, lowerEdge, rightEdge, upperEdge);
            postfit_error_boxes[i_bin - 1].SetFillStyle(3004);
            postfit_error_boxes[i_bin - 1].SetFillColor(kGray + 3);
            postfit_error_boxes[i_bin - 1].Draw("same");
        }


        cv->Draw();
        cv->SaveAs("test.png");
    }
}

void hf003_example()
{
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
    meas.SetLumiRelErr( 0.00 );
    meas.SetExportOnly( false );
    meas.SetBinHigh( 2 );

    // Create a channel

    Channel chan( "channel1" );
    chan.SetData( "data", InputFile );
    chan.SetStatErrorConfig( 0.05, "Poisson" );

    // Now, create some samples


    // Create the signal sample
    Sample signal( "signal", "signal", InputFile );
    signal.ActivateStatError();
    signal.AddNormFactor("SigXsecOverSM", 1, 0, 3 );
    signal.AddOverallSys("syst1", 0.95, 1.05);
    chan.AddSample( signal );

    // Background 1
    Sample background1("background1", "background1", InputFile );
    // Input file name can be skipped, in such case same input file as for nominal histogram will be used
    background1.AddHistoSys("weighted_modeling",
                                "weighted_modeling_up", InputFile, "",
                                "weighted_modeling_down", InputFile, "");
    background1.AddHistoSys("modeling",
                                "modeling_up", InputFile, "", 
                                "modeling_down", InputFile, "");
    background1.AddOverallSys("syst1", 0.95, 1.05);
    background1.ActivateStatError();
    chan.AddSample( background1 );



    // Done with this channel
    // Add it to the measurement:
    meas.AddChannel( chan );

    // Collect the histograms from their files,
    meas.CollectHistograms();

    // Simple histogram can be drawn using only measurement
    DrawHistogramUsingMeasurement(meas); 

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

    DrawHistogramUsingRooWorkspace(ws.get(), meas, result.get());
    // Draw prefit and postfit histograms including uncertainties calculation

}
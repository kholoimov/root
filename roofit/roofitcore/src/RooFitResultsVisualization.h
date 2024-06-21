#include "TCanvas.h"
#include "TString.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include <vector>
#include "TBox.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLine.h"
#include "TStyle.h"



class RooFitResultsVisualization
{
public:
    RooFitResultsVisualization(RooWorkspace* ws) : fWorkspace(ws) {}
    virtual ~RooFitResultsVisualization() {}

    

private:
    RooWorkspace* fWorkspace = nullptr;
    mutable std::vector<TString> names;
    mutable std::vector<double> values;
    mutable std::vector<double> errors;

};

void makeData()
{
  TFile* example = new TFile("hf002_input_data.root","RECREATE");
  TH1F* data = new TH1F("data","data", 4, 200, 600);
  TH1F* signal = new TH1F("signal","signal histogram (pb)", 4, 200, 600);
  TH1F* background1 = new TH1F("background1","background 1 histogram (pb)", 4, 200, 600);
//   TH1F* background2 = new TH1F("background2","background 2 histogram (pb)", 4, 200, 600);
//   TH1F* statUncert = new TH1F("background1_statUncert", "statUncert", 4, 200, 600);

  TH1F* weighted_modeling_up = new TH1F("weighted_modeling_up","weighted_modeling_up histogram (pb)", 4, 200, 600);
  TH1F* weighted_modeling_down = new TH1F("weighted_modeling_down","weighted_modeling_down histogram (pb)", 4, 200, 600);
  TH1F* modeling_up = new TH1F("modeling_up","modeling_up histogram (pb)", 4, 200, 600);
  TH1F* modeling_down = new TH1F("modeling_down","modeling_down histogram (pb)", 4, 200, 600);
  

  // run with 1 pb
  data->SetBinContent(1,112);
  data->SetBinContent(2,112);
  data->SetBinContent(3,120);
  data->SetBinContent(4,65);

  signal->SetBinContent(1,0);
  signal->SetBinContent(2,2);
  signal->SetBinContent(3,21);
  signal->SetBinContent(4,24);

  background1->SetBinContent(1,112);
  background1->SetBinContent(2,130);
  background1->SetBinContent(3,85);
  background1->SetBinContent(4,55);


  weighted_modeling_up->SetBinContent(1,98);
  weighted_modeling_up->SetBinContent(2,152);
  weighted_modeling_up->SetBinContent(3,132);
  weighted_modeling_up->SetBinContent(4,100);

  weighted_modeling_down->SetBinContent(1,79);
  weighted_modeling_down->SetBinContent(2,89);
  weighted_modeling_down->SetBinContent(3,62);
  weighted_modeling_down->SetBinContent(4,39);

  modeling_up->SetBinContent(1,54);
  modeling_up->SetBinContent(2,85);
  modeling_up->SetBinContent(3,90);
  modeling_up->SetBinContent(4,82);

  modeling_down->SetBinContent(1,170);
  modeling_down->SetBinContent(2,171);
  modeling_down->SetBinContent(3,85);
  modeling_down->SetBinContent(4,30);
  // A small statistical uncertainty

  example->Write();
  example->Close();
  ////////////////////// 
}

void makeExample002()
{
    makeData();
}
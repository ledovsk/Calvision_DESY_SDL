#include <TROOT.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <Riostream.h>
#include <TVirtualFFT.h>
#include <TLegend.h>

constexpr unsigned int N_Channels = 2 * 10;
constexpr unsigned int N_Samples = 1024;

TGraph *grCerenkov;
TGraph *grScint;

pair<double,double> baseline(TGraph*, double, double);
double funcSDL(double*, double*);
TGraph *grGetScintillation(TGraph*, double, double, double, double, double);







class TreeReader {
public:
  TreeReader(TString filename, TString treename="tree") {
    tree_ = new TChain(treename);
    tree_->Add(filename);

    /*
    std::cout << "The following files were added to the TChain:\n";
    const auto& fileElements = *tree_->GetListOfFiles();
    for (TObject const* op : fileElements) {
      auto chainElement = static_cast<const TChainElement*>(op);
      std::cout << chainElement->GetTitle() << "\n";
    }
    
    tree_->Print();
    */
    tree_->SetBranchAddress("horizontal_interval", &horizontal_interval_);
    tree_->SetBranchAddress("channels", channels_[0].data());
    tree_->SetBranchAddress("trigger", trigger_[0].data());

  }
  
  UInt_t num_entries() const {
    //return 5000; // std::min(5000, tree_->GetEntries());
    return tree_->GetEntries();
  }

  void get_entry(UInt_t i) {
    tree_->GetEntry(i);
  }

  const std::array<Double_t, N_Samples> time() const {
    std::array<Double_t, N_Samples> t;
    for (unsigned int i = 0; i < N_Samples; i++) {
      t[i] = static_cast<Double_t>(i) * horizontal_interval_;
    }
    return t;
  }

  std::array<Double_t, N_Samples> voltages(UInt_t channel) const {
    std::array<Double_t, N_Samples> volts;
    for (unsigned int i = 0; i < N_Samples; i++) {
      volts[i] = static_cast<Double_t>(channels_[channel][i]);
    }
    return volts;
  }

  std::array<Double_t, N_Samples> trigger(UInt_t channel) const {
    std::array<Double_t, N_Samples> volts;
    for (unsigned int i = 0; i < N_Samples; i++) {
      volts[i] = static_cast<Double_t>(trigger_[channel][i]);
    }
    return volts;
  }

  
  Float_t horizontal_interval() const { return horizontal_interval_; }

  UInt_t get_bin(Double_t t) const {
    return static_cast<UInt_t>(t / horizontal_interval_);
  }

private:
  TChain* tree_;

  Float_t horizontal_interval_;
  std::array<std::array<Float_t, N_Samples>, N_Channels> channels_;
  std::array<std::array<Float_t, N_Samples>, 2> trigger_;
  std::array<Double_t, N_Samples> times_;
};



class RecoTRG{
public:
    RecoTRG(TGraph *gr){

      pair<double,double> p = baseline(gr,   0, 100);
      pair<double,double> a = baseline(gr, 400, 1000);
      ped_ = p.first;
      rms_ = p.second;
      amp_ = a.first - p.first;

      int     N = gr->GetN();
      double *X = gr->GetX();
      double *Y = gr->GetY();
      // Timing of the threshold using linear interpolation
      double threshold = ped_ + 0.5 * amp_;
      time_ = 0;
      for(int i=1; i<N; i++){
          if( Y[i-1] <= threshold && Y[i] > threshold){
              time_ = X[i-1] + (threshold - Y[i-1])/(Y[i] - Y[i-1]) * (X[i] - X[i-1]);
              break;
          }
      }
    }
            
    double ped()  { return ped_; }
    double rms()  { return rms_; }
    double amp()  { return amp_; }
    double time() { return time_; }
        

private:
    double ped_;
    double rms_;
    double amp_;
    double time_;
};



pair<double,double> baseline(TGraph *gr, double tmin=0, double tmax=1000.)
{
  int     N = gr->GetN();
  double *X = gr->GetX();
  double *Y = gr->GetY();
  
  std::vector<double> v;
  for(int i=0; i<N; i++){
    if(X[i] > tmin && X[i] < tmax)
      v.push_back(Y[i]);
  }
  
  double vmin = 0.;
  double vmax = 0.;
  std::sort(v.begin(),v.end());
  unsigned int nInInterval = int(0.68*v.size());
  double minInterval = 1e+9;
  for (unsigned int i=0; i<v.size(); i++) {
    double interval_size=minInterval;
    if ((i+nInInterval)<v.size()) interval_size = (v[i+nInInterval]-v[i]);
    if (interval_size < minInterval){
      minInterval = interval_size;
      vmin = v[i];
      vmax = vmin + minInterval;
    }
  }
  
  pair<double,double> result;
  result.first = 0.5 * (vmax + vmin);
  result.second = 0.5 * (vmax - vmin);
  return result;
}




double funcSDL(double *x, double *par)
{
    double t = x[0] - par[0];
    double f = 0;
    if(t > -70){
        f += par[1] * grCerenkov->Eval(t) + par[2] * grScint->Eval(t);
    }
    return f;
}



TGraph *grGetScintillation(TGraph *gr0, double tauS=300, double tauM=300, double tauL=300, double fracS=0.0, double fracM=0.0)
{
    /*
     * Build average 1pe pulse shape for scintillation photo-electron based on average 1pe pulse shape for Cerenkov one
     * Assume that scintillation has three decay components: short, medium and long with decay times tauS, tauM and tauL, respectively
     * Fractions of photo-electrons with tauS and tauM are fracS and fracM, respectively
     */
    
    int     N = gr0->GetN();
    double *X = gr0->GetX();
    double *Y = gr0->GetY();
    
    TF1 *fDecay = new TF1("fDecay","[3]*exp(-x/[0])/[0]+[4]*exp(-x/[1])/[1]+(1-[3]-[4])*exp(-x/[2])/[2]",0,2000);
    fDecay->SetParameters(tauS,tauM,tauL,fracS,fracM);
    std::vector<double> Z;
    for(int i=0; i<N; i++) Z.push_back(0);
    
    double dt = 0.1;
    double tNow = 0.5 * dt;
    while(tNow < 1000){
        double w = dt * fDecay->Eval(tNow);
        for(int i=0; i<N; i++){
            double t = X[i] - tNow;
            if(t>-70){
                Z[i] += w * gr0->Eval(t);
            }
        }
        tNow += dt;
    }
    TGraph *gr = new TGraph();
    for(int i=0; i<N; i++){
        gr->SetPoint(i, X[i], Z[i]);
    }
    delete fDecay;
    return gr;
}



void plot_DRS_channel(int event=0, int ch=0, TString filename="../data/run_374/outfile_LG.root")
{
  TreeReader tree(filename);
  tree.get_entry(event);
  
  TGraph *grWF = new TGraph(N_Samples, tree.time().data(), tree.voltages(ch).data());
  grWF->GetXaxis()->SetTitle("t (ns)");
  grWF->GetYaxis()->SetTitle("a (mV)");
  grWF->Draw("APL");
}



void plot_DRS_trigger(int event=0, TString filename="../data/run_374/outfile_LG.root")
{
  TreeReader tree(filename);
  tree.get_entry(event);
  
  TGraph *grWF = new TGraph(N_Samples, tree.time().data(), tree.trigger(0).data());
  grWF->GetXaxis()->SetTitle("t (ns)");
  grWF->GetYaxis()->SetTitle("a (mV)");
  grWF->Draw("APL");
  
  RecoTRG *trg  = new RecoTRG(grWF);
  cout << " Baseline RMS = " << trg->rms() << endl;
  cout << " Amplitude    = " << trg->amp() << endl;
  cout << " Time @ 50%   = " << trg->time() << endl;
  
}



void reco_WF_SDL(int ievt=0, int ch=0)
{
    /*
     * Example of reconstruction for DSB crystal:
     * Run = 200
     * Channels 0, 1, 2, and 3 are in the front
     * Channels 4 and 7 are in the rear
     * Channels 5 and 6 are discarded from the analysis (dead)
     */
    if(!(ch==0 || ch==1 || ch==2 || ch==3 || ch==4 || ch==7)) return;
    TString filename = "../data/run_200/outfile_LG.root";

    
    // SPR for Cerenkov
    TFile *fSPR = new TFile("SPR_SDL.root");
    grCerenkov = (TGraph*)fSPR->Get(Form("grSPR_SDL_ch%d_rear",ch));

    // SPR for Scintillation
    // DSB has two decay components 100ns (13%) and 500ns (87%)
    grScint    = grGetScintillation(grCerenkov, 30,100,500,0.0,0.13);

    // Function for component fit
    TF1 *fFit = new TF1("fFit",funcSDL,-100,800,3);
    fFit->SetLineColor(kOrange+10);
    fFit->SetParName(0,"time jitter");
    fFit->SetParName(1,"Nc x A1pe");
    fFit->SetParName(2,"Ns x A1pe");
  
  
    TreeReader tree(filename);
    tree.get_entry(ievt);

    // Find trigger timing
    TGraph *grTr = new TGraph(N_Samples, tree.time().data(), tree.trigger(0).data());
    RecoTRG *trig = new RecoTRG( grTr );
    delete grTr;
    
    // DRS waveform
    TGraph *grWF = new TGraph(N_Samples, tree.time().data(), tree.voltages(ch).data());

    
    // SDL waveform with 2ns delay, adjusted for trigger timing
    int     N = grWF->GetN();
    double *X = grWF->GetX();
    double *Y = grWF->GetY();
    TGraph *grSDL = new TGraph();
    for(int it=2; it<N; it++){
            grSDL->SetPoint(it-2, X[it] - trig->time(), Y[it-2] - Y[it]);
    }
    fFit->SetParameters(0, 10, 500.);
    fFit->SetParLimits(0, -2, 2);
    fFit->SetParLimits(1, 0, 1000);
    fFit->SetParLimits(2, 0, 100000.);
    grSDL->Fit("fFit","WN","",-60,50);
    
    grSDL->GetXaxis()->SetLimits(-100,50);
    grSDL->GetXaxis()->SetTitle("t - t_{ TRG} (ns)");
    grSDL->GetYaxis()->SetTitle("a(t) - a(t-2ns) (mV)");
          
    TF1 *fC = new TF1("fC",funcSDL,-100,800,3);
    TF1 *fS = new TF1("fS",funcSDL,-100,800,3);
    fC->SetLineColor(kAzure+6);
    fS->SetLineColor(kSpring-1);

    grSDL->Draw("APL");
    fC->SetParameters(fFit->GetParameter(0), fFit->GetParameter(1), 0);
    fS->SetParameters(fFit->GetParameter(0), 0, fFit->GetParameter(2));
    fFit->SetNpx(10000);
    fC->SetNpx(10000);
    fS->SetNpx(10000);
    fFit->SetLineWidth(1);
    fC->SetLineWidth(1);
    fS->SetLineWidth(1);
    
    fFit->Draw("same");
    fC->Draw("same");
    fS->Draw("same");
         
    TLegend *leg = new TLegend(0.54,0.84,0.97,0.99,NULL,"brNDC");
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.033);
    leg->AddEntry(grSDL,"Measurements","pl");
    leg->AddEntry(fFit,"Fit = Cerenkov + Scintillation","l");
    leg->AddEntry(fC,"Cerenkov","l");
    leg->AddEntry(fS,"Scintillation","l");
    leg->Draw(); 
    
}


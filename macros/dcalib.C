enum type_t { alcor_hit = 1, trigger_tag = 9, start_spill = 7, end_spill = 15 };

struct hit_t {
  int fifo, type, counter, column, pixel, tdc, rollover, coarse, fine;
  double time;
};

void
dcalib(const std::string infilename, const std::string outfilename, int channel)
{
  
  /** open input file, get tree and
      connect variables to branches **/
  auto fin = TFile::Open(infilename.c_str());
  auto tin = (TTree *)fin->Get("alcor");
  auto nev = tin->GetEntries();
  hit_t hit;
  tin->SetBranchAddress("fifo", &hit.fifo);
  tin->SetBranchAddress("type", &hit.type);
  tin->SetBranchAddress("counter", &hit.counter);
  tin->SetBranchAddress("column", &hit.column);
  tin->SetBranchAddress("pixel", &hit.pixel);
  tin->SetBranchAddress("tdc", &hit.tdc);
  tin->SetBranchAddress("rollover", &hit.rollover);
  tin->SetBranchAddress("coarse", &hit.coarse);
  tin->SetBranchAddress("fine", &hit.fine);
  
  /** minimisation function **/

  auto chi2 = [&](const double *par, TH1 *hDelta = nullptr) {
    double f = 0.;

    bool has_first_hit_in_spill = false;
    hit_t phit;
    int n = 0;
    for (int iev = 0; iev < nev; ++iev) {
      tin->GetEntry(iev);

      /** reset when start of spill is detected **/
      if (hit.type == start_spill) {
	has_first_hit_in_spill = false;
	continue;
      }

      /** must be ALCOR hit **/
      if (hit.type != alcor_hit)
	continue;

      /** from the desired channel **/
      if (hit.column * 4 + hit.pixel != channel)
	continue;

      /** calculate corrected time **/
      int coarse = hit.coarse + 32768 * hit.rollover;
      double corr = par[hit.tdc] + (double)hit.fine * par[hit.tdc + 4];
      hit.time = (double)coarse - corr;

      /** is this is the first hit in the spill, store it and continue **/
      if (!has_first_hit_in_spill) {
	phit = hit;
	has_first_hit_in_spill = true;
	continue;
      }

      /** compute time difference wrt. previous hit and wrt. period **/
      double delta = hit.time - phit.time - par[8];

      /** this is a safety measure, because sometimes there is some crap **/
      if ( fabs(coarse - (phit.coarse + 32768 * phit.rollover) - 320) < 10. ) {
	f += delta * delta;
	++n;
      }
      
      if (hDelta) hDelta->Fill(delta);

      phit = hit;    
    }
    std::cout << " --- minimising: sigma = " << std::sqrt(f / n) / sqrt(2.) * 3.125 << " ns " << std::endl;
    return f;
  };

  ROOT::Math::Functor fcn(chi2, 9);
  ROOT::Fit::Fitter fitter;
  double pStart[9] = { 0.5, 0.5, 0.5, 0.5, 0.0156, 0.0156, 0.0156, 0.0156, 320. };
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("off_0");
  fitter.Config().ParSettings(1).SetName("off_1");
  fitter.Config().ParSettings(2).SetName("off_2");
  fitter.Config().ParSettings(3).SetName("off_3");
  fitter.Config().ParSettings(4).SetName("iif_0");
  fitter.Config().ParSettings(5).SetName("iif_1");
  fitter.Config().ParSettings(6).SetName("iif_2");
  fitter.Config().ParSettings(7).SetName("iif_3");
  fitter.Config().ParSettings(8).SetName("period");

  fitter.Config().ParSettings(0).Fix();

  bool ok = fitter.FitFCN();
  auto result = fitter.Result();
  result.Print(std::cout);
  double par[9], pare[9];
  auto hParam = new TH1F("hParam", "", 9, 0., 9.);
  for (int ipar = 0; ipar < 9; ++ipar) {
    par[ipar] = result.Parameter(ipar);
    pare[ipar] = result.ParError(ipar);
    hParam->SetBinContent(ipar + 1, par[ipar]);
    hParam->SetBinError(ipar + 1, pare[ipar]);
  }

  /** draw **/

  auto hDelta = new TH1F("hDelta", "", 50, -0.5, 0.5);
  auto res = chi2(par, hDelta);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  hParam->Write();
  hDelta->Write();
  fout->Close();
  
}

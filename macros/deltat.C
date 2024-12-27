/**
   This macro measures the time distance between two consecutive ALCOR hits in the same channel.
   It also allows the user to upload ALCOR TDC fine calibration parameters of the channel.
   
   The fine correction of each of the four TDCs on one ALCOR channel is defined as
   
   corr = offset_i + slope_i * fine
   
   where i is the TDC index.
   The calibration parameters are stores in a double[8] array as follows
   
   { offset_0, offset_1, offset_2, offset_3, slope_0, slope_1, slope_2, slope_3 }
    
 **/

enum type_t {
  alcor_hit = 1, trigger_tag = 9, start_spill = 7, end_spill = 15
};

struct hit_t {
  int fifo, type, counter, column, pixel, tdc, rollover, coarse, fine;
  double time;
};

void
deltat(const std::string filename, int channel, const std::string calibfilename = "")
{
  
  /** open input file and get tree **/
  std::cout << " --- opening input file: " << filename << std::endl;
  auto fin = TFile::Open(filename.c_str());
  if (!fin || !fin->IsOpen()) {
    std::cerr << " --- could not open input file " << std::endl;
    return;
  }
  auto tin = (TTree *)fin->Get("alcor");
  if (!tin) {
    std::cerr << " --- could not find \'alcor\' tree in input file " << std::endl;
    return;
  }
  auto nev = tin->GetEntries();
  std::cout << " --- found " << nev << " entries in the input tree " << std::endl;

  /** create hit and connect variables to branches **/
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

  /** open calibration file and load calibration parameters **/
  double par[8] = {0.};
  if (!calibfilename.empty()) {
    std::cout << " --- loading calibration: " << calibfilename << std::endl;
    auto fcal = TFile::Open(calibfilename.c_str());
    auto hParam = (TH1 *)fcal->Get("hParam");
    for (int ipar = 0; ipar < 8; ++ipar)
      par[ipar] = hParam->GetBinContent(ipar + 1);
  }
  
  auto hDelta = new TH1F("hDelta", "", 50000, 0., 1000.);
  
  /** loop over tree entries **/
  std::cout << " --- loop over tree entries " << std::endl;
  bool has_first_hit_in_spill = false;
  hit_t phit;
  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);

    /** reset when start of spill is detected **/
    if (hit.type == start_spill) {
      std::cout << " --- start of spill " << std::endl;
      has_first_hit_in_spill = false;
      continue;
    }

    /** must be a ALCOR hit **/
    if (hit.type != alcor_hit)
      continue;

    /** from the desired channel **/
    if (hit.column * 4 + hit.pixel != channel)
      continue;

    /** calculate corrected time **/
    double coarse = (double)hit.coarse + 32768. * (double)hit.rollover;
    double fine = par[hit.tdc] + (double)hit.fine * par[hit.tdc + 4];
    hit.time = coarse - fine;

    /** is this is the first hit in the spill, store it and continue **/
    if (!has_first_hit_in_spill) {
      phit = hit;
      has_first_hit_in_spill = true;
      continue;
    }

    /** compute time difference wrt. previous hit **/
    auto delta = hit.time - phit.time;
    hDelta->Fill(delta);

    /** store current hit **/
    phit = hit;
  }

  hDelta->Draw();
}

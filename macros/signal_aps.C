enum type_t { alcor_hit = 1, trigger_tag = 9, start_spill = 7, end_spill = 15 };

typedef std::array<int, 2> frame_t;
const int frame_size = 256;

const bool afterpulse_suppression = false;
const double dead_time = 32.;

struct hit_t {
  int fifo, type, counter, column, pixel, tdc, rollover, coarse, fine;
  double time;
  frame_t frame;
  static bool compare(hit_t &a, hit_t &b) { return a.time < b.time; };
};

typedef std::map<frame_t, std::vector<hit_t>> hitmap_t;

void collect_hits(hitmap_t &hits, const std::string filename, int channel, const std::string calibfilename, bool aps = false);

void
signal_aps(const std::string reference_filename, int reference_channel, const std::string reference_calibfilename,
	   const std::string target_filename, int target_channel, const std::string target_calibfilename,
	   const std::string outfilename)
{

  /** collect hits from reference and target data **/
  hitmap_t reference_hits, target_hits;
  collect_hits(reference_hits, reference_filename, reference_channel, reference_calibfilename);
  collect_hits(target_hits, target_filename, target_channel, target_calibfilename, afterpulse_suppression);

  auto hDelta = new TH1F("hDelta", ";hit - reference time (clock cycles);normalised counts", 10000, -100, 100);
  int ntriggers = 0;

  /** loop over reference frames **/
  for (auto &[frame, rhits] : reference_hits) {

    /** loop over reference hits in frame **/
    for (auto rhit : rhits) {
      ++ntriggers;

      /** loop over target hits in frame **/
      for (auto thit : target_hits[frame]) {
	double delta = thit.time - rhit.time;
	hDelta->Fill(delta);
      }
      
    }
  }

  /** normalise to number of triggers **/
  hDelta->Sumw2();
  hDelta->Scale(1. / ntriggers);

  /** apply geometrical acceptance function to obtain corrected distribution **/
  TF1 *facc = new TF1("facc", "([0] + x) * (x >= -[0] && x < 0) / [0] + ([0] - x) * (x >= 0 && x <= [0]) / [0]", -100, 100);
  facc->SetParameter(0, frame_size);
  auto hDelta_corr = (TH1 *)hDelta->Clone("hDelta_corr");
  hDelta_corr->Divide(facc);

  std::cout << " --- writing output: " << outfilename << std::endl;
  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  hDelta->Write();
  hDelta_corr->Write();
  fout->Close();

}

void
collect_hits(hitmap_t &hits, const std::string filename, int channel, const std::string calibfilename, bool aps)
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
    fcal->Close();
  }

  /** temporary hit vector **/
  std::vector<hit_t> _hits;
  
  /** loop over tree entries **/
  std::cout << " --- loop over tree entries " << std::endl;
  int spill = 0;
  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);

    /** end of spill **/
    if (hit.type == end_spill) {
      std::cout << " --- end of spill " << std::endl;

      /** sort hits in increasing time order **/
      std::sort(_hits.begin(), _hits.end(), hit_t::compare);
      
      /** push hits to the output hit map **/
      hit_t phit;
      bool has_phit = false;
      for (auto &_hit : _hits) {
	if (aps && has_phit && (_hit.time - phit.time) < dead_time)
	  continue;
	hits[_hit.frame].push_back(_hit);
	phit = _hit;
	has_phit = true;
      }
      
      ++spill;
      continue;
    }
    
    /** reset start of spill **/
    if (hit.type == start_spill) {
      std::cout << " --- start of spill " << std::endl;
      _hits.clear();
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

    /** define frame **/
    int subframe = coarse / frame_size;
    frame_t frame = {spill, subframe};
    hit.frame = frame;

    /** push hit **/
    _hits.push_back(hit);

  }
  std::cout << " --- collected " << hits.size() << " frames " << std::endl;
  
  /** close input file **/
  fin->Close();

}


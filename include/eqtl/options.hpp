/* eqtl — options / CLI */
#pragma once
#include <string>
#include <vector>

namespace eqtl {

enum class Mode { Cis, Trans, All, Gw };
enum class Model { Lm, Glm, Lmm, Glmm };
enum class MissHand { Filter, Impute };

struct Options {
  // Genotypes: exactly one of vcf or bfile (PLINK prefix → .bed/.bim/.fam)
  std::string vcf;
  std::string bfile;
  std::string pheno;
  std::string out = "eqtl_out";
  std::string gff;
  std::string covar;
  std::string grm;
  std::string gff_id_key; // empty = try ID,Name,gene_id

  bool make_grm = false;
  Mode mode = Mode::All;
  std::vector<Model> models{Model::Lmm};

  int window = 1000000;
  double pval_cis = 1e-5;
  double pval_trans = 1e-5;

  MissHand miss = MissHand::Impute;
  // drop SNP if missing fraction among analysis samples > max_miss
  // filter default 0 => any missing drops; impute still respects max_miss before fill
  double max_miss = 0.8;
  // keep if maf_min <= MAF <= 1-maf_min on analysis samples (non-missing); 0 = off
  double maf = 0.0;
  bool fast = false;
  int threads = 1;

  int perm = 0; // gene-level permutations; 0 = off
  int seed = -1; // -1 = unset
  bool disable_beta_approx = false;

  bool help = false;
  bool version = false;

  // fission subcommand
  bool   run_fission    = false;
  int    peer_factors   = 10;
  double fission_epsilon = 0.5;
  int    fission_max_iter = 1000;
  double fission_tol    = 1e-3;

  bool use_bfile() const { return !bfile.empty(); }
};

// returns 0 ok, 1 error (message on stderr)
int parse_options(int argc, char** argv, Options& opt);
void print_help();
void print_version();

std::string mode_str(Mode m);
std::string model_str(Model m);
const char* miss_str(MissHand m);

} // namespace eqtl

#include "eqtl/options.hpp"
#include "eqtl/util.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <getopt.h>

namespace eqtl {

static Model parse_one_model(const std::string& s) {
  if (s == "lm") return Model::Lm;
  if (s == "glm") return Model::Glm;
  if (s == "lmm") return Model::Lmm;
  if (s == "glmm") return Model::Glmm;
  die("invalid --model value: " + s + " (lm|glm|lmm|glmm)");
  return Model::Lmm;
}

static std::vector<Model> parse_models(const std::string& s) {
  std::vector<Model> out;
  for (auto& p : split_char(s, ',')) {
    auto t = trim(p);
    if (!t.empty()) out.push_back(parse_one_model(t));
  }
  if (out.empty()) die("empty --model");
  return out;
}

static Mode parse_mode(const std::string& s) {
  if (s == "cis") return Mode::Cis;
  if (s == "trans") return Mode::Trans;
  if (s == "all") return Mode::All;
  if (s == "gw" || s == "gwas") return Mode::Gw;
  die("invalid --mode: " + s + " (cis|trans|all|gw)");
  return Mode::All;
}

std::string mode_str(Mode m) {
  switch (m) {
    case Mode::Cis: return "cis";
    case Mode::Trans: return "trans";
    case Mode::All: return "all";
    case Mode::Gw: return "gw";
  }
  return "all";
}

std::string model_str(Model m) {
  switch (m) {
    case Model::Lm: return "lm";
    case Model::Glm: return "glm";
    case Model::Lmm: return "lmm";
    case Model::Glmm: return "glmm";
  }
  return "lmm";
}

const char* miss_str(MissHand m) {
  return m == MissHand::Filter ? "filter" : "impute";
}

void print_version() {
  std::cout << "eqtl 0.1.0\n";
}

void print_help() {
  std::cout
    << "Program: eqtl (expression QTL mapping)\n"
    << "Version: 0.1.0\n"
    << "\n"
    << "Usage:   eqtl [options]\n"
    << "\n"
    << "Input options:\n"
    << "    -v, --vcf FILE         genotypes VCF/BCF (GT)\n"
    << "    -e, --pheno FILE       phenotype matrix (col1=sample)\n"
    << "    -g, --gff FILE         GFF3 gene annotation (TSS)\n"
    << "    -c, --covar FILE       covariates\n"
    << "    -k, --grm PREFIX       relatedness matrix (.grm.id/.grm.bin)\n"
    << "        --gff-id-key STR   GFF attribute for gene id\n"
    << "\n"
    << "Analysis options:\n"
    << "        --make-grm         write relatedness matrix and exit\n"
    << "    -m, --mode STR         cis|trans|all|gw  [all]\n"
    << "        --model STR        lm|glm|lmm|glmm[,...]  [lmm]\n"
    << "    -w, --window INT       cis window around TSS (bp)  [1000000]\n"
    << "        --pval-cis FLOAT   cis output p threshold  [1e-5]\n"
    << "        --pval-trans FLOAT trans/gw output p threshold  [1e-5]\n"
    << "        --miss-hand STR    filter|impute missing GT  [filter]\n"
    << "        --max-miss FLOAT   drop SNP if missing fraction > value  [0]\n"
    << "        --fast             share variance/dispersion params per gene\n"
    << "\n"
    << "Permutation options:\n"
    << "        --perm INT         gene-level permutations  [10000]\n"
    << "        --seed INT         permutation seed\n"
    << "        --disable-beta-approx  omit beta-approximated p\n"
    << "\n"
    << "Output options:\n"
    << "    -o, --out PREFIX       output prefix  [eqtl_out]\n"
    << "    -t, --thread INT       threads  [1]\n"
    << "    -h, --help             help\n"
    << "        --version          version\n"
    << "\n"
    << "Examples:\n"
    << "    eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff --model lm --mode cis --perm 0 -o out\n"
    << "    eqtl -v data/smoke.vcf.gz --make-grm -o data/smoke\n"
    << "    eqtl -v in.vcf.gz -e pheno.tsv -g genes.gff --model lmm,glm --fast -o run1\n";
}

int parse_options(int argc, char** argv, Options& opt) {
  if (argc == 1) {
    print_help();
    return 2;
  }

  static struct option long_opts[] = {
      {"vcf", required_argument, 0, 'v'},
      {"pheno", required_argument, 0, 'e'},
      {"out", required_argument, 0, 'o'},
      {"gff", required_argument, 0, 'g'},
      {"covar", required_argument, 0, 'c'},
      {"grm", required_argument, 0, 'k'},
      {"gff-id-key", required_argument, 0, 1001},
      {"make-grm", no_argument, 0, 1002},
      {"mode", required_argument, 0, 'm'},
      {"model", required_argument, 0, 1003},
      {"window", required_argument, 0, 'w'},
      {"pval-cis", required_argument, 0, 1004},
      {"pval-trans", required_argument, 0, 1005},
      {"miss-hand", required_argument, 0, 1006},
      {"max-miss", required_argument, 0, 1012},
      {"fast", no_argument, 0, 1007},
      {"thread", required_argument, 0, 't'},
      {"threads", required_argument, 0, 't'},
      {"perm", required_argument, 0, 1008},
      {"permutations", required_argument, 0, 1008},
      {"seed", required_argument, 0, 1009},
      {"disable-beta-approx", no_argument, 0, 1010},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 1011},
      {0, 0, 0, 0}};

  int c;
  int idx;
  while ((c = getopt_long(argc, argv, "v:e:o:g:c:k:m:w:t:h", long_opts, &idx)) != -1) {
    switch (c) {
      case 'v': opt.vcf = optarg; break;
      case 'e': opt.pheno = optarg; break;
      case 'o': opt.out = optarg; break;
      case 'g': opt.gff = optarg; break;
      case 'c': opt.covar = optarg; break;
      case 'k': opt.grm = optarg; break;
      case 'm': opt.mode = parse_mode(optarg); break;
      case 'w': opt.window = std::atoi(optarg); break;
      case 't': opt.threads = std::max(1, std::atoi(optarg)); break;
      case 'h': opt.help = true; break;
      case 1001: opt.gff_id_key = optarg; break;
      case 1002: opt.make_grm = true; break;
      case 1003: opt.models = parse_models(optarg); break;
      case 1004: opt.pval_cis = std::atof(optarg); break;
      case 1005: opt.pval_trans = std::atof(optarg); break;
      case 1006:
        if (std::string(optarg) == "filter") opt.miss = MissHand::Filter;
        else if (std::string(optarg) == "impute") opt.miss = MissHand::Impute;
        else die("invalid --miss-hand (filter|impute)");
        break;
      case 1007: opt.fast = true; break;
      case 1008: opt.perm = std::atoi(optarg); break;
      case 1009: opt.seed = std::atoi(optarg); break;
      case 1010: opt.disable_beta_approx = true; break;
      case 1011: opt.version = true; break;
      case 1012:
        opt.max_miss = std::atof(optarg);
        if (opt.max_miss < 0.0 || opt.max_miss > 1.0) die("--max-miss must be in [0,1]");
        break;
      default:
        return 1;
    }
  }

  if (opt.help) {
    print_help();
    return 2;
  }
  if (opt.version) {
    print_version();
    return 2;
  }

  if (opt.vcf.empty()) {
    std::cerr << "[E] missing -v/--vcf\n";
    return 1;
  }
  if (!opt.make_grm && opt.pheno.empty()) {
    std::cerr << "[E] missing -e/--pheno (or use --make-grm)\n";
    return 1;
  }
  if (opt.window < 0) die("--window must be >= 0");
  return 0;
}

} // namespace eqtl

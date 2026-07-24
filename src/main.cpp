#include "eqtl/options.hpp"
#include "eqtl/fission.hpp"
#include "eqtl/scan.hpp"
#include "eqtl/util.hpp"

int main(int argc, char** argv) {
  eqtl::Options opt;
  const int pr = eqtl::parse_options(argc, argv, opt);
  if (pr == 2) {
    return 0;  // help/version
  }
  if (pr != 0) {
    return 1;
  }
  try {
    if (opt.run_fission) return eqtl::run_fission(opt);
    if (opt.make_grm) return eqtl::run_make_grm(opt);
    return eqtl::run_eqtl(opt);
  } catch (const std::exception& ex) {
    eqtl::die(std::string("exception: ") + ex.what());
  }
  return 0;
}

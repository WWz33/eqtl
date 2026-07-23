/* eqtl — scan + write */
#pragma once
#include "eqtl/options.hpp"
#include "eqtl/pheno.hpp"
#include "eqtl/annot.hpp"
#include "eqtl/grm.hpp"
#include "eqtl/vcf_session.hpp"
#include "eqtl/models.hpp"

namespace eqtl {

int run_make_grm(const Options& opt);
int run_eqtl(const Options& opt);

} // namespace eqtl

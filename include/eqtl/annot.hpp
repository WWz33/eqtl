/* eqtl — GFF gene TSS (via third_party/gffsub) */
#pragma once
#include <string>
#include <vector>
#include <unordered_map>

namespace eqtl {

struct GeneLoc {
  std::string id;
  std::string chrom;
  int64_t tss = 0;
  char strand = '+';
  bool ok = false;
};

// Parse GFF3 genes. id_key empty: try ID then Name then gene_id.
// multi-hit: warn skip; 0-hit genes remain missing.
std::unordered_map<std::string, GeneLoc> load_gff_tss(
    const std::string& path, const std::string& id_key);

} // namespace eqtl

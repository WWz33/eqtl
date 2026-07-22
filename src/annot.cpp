#include "eqtl/annot.hpp"
#include "eqtl/util.hpp"
#include <fstream>
#include <unordered_set>

namespace eqtl {

static std::string attr_value(const std::string& attrs, const std::string& key) {
  // GFF3: key=value;key2=value2
  size_t pos = 0;
  while (pos < attrs.size()) {
    size_t eq = attrs.find('=', pos);
    if (eq == std::string::npos) break;
    std::string k = attrs.substr(pos, eq - pos);
    size_t sc = attrs.find(';', eq + 1);
    std::string v = (sc == std::string::npos) ? attrs.substr(eq + 1)
                                              : attrs.substr(eq + 1, sc - eq - 1);
    if (k == key) return v;
    if (sc == std::string::npos) break;
    pos = sc + 1;
  }
  return {};
}

std::unordered_map<std::string, GeneLoc> load_gff_tss(
    const std::string& path, const std::string& id_key) {
  std::ifstream in(path);
  if (!in) die("cannot open gff: " + path);
  std::unordered_map<std::string, GeneLoc> out;
  std::unordered_set<std::string> multi;
  std::string line;
  int n_gene = 0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    auto t = split_char(line, '\t');
    if (t.size() < 9) continue;
    if (t[2] != "gene") continue;
    ++n_gene;
    std::string id;
    if (!id_key.empty()) {
      id = attr_value(t[8], id_key);
    } else {
      id = attr_value(t[8], "ID");
      if (id.empty()) id = attr_value(t[8], "Name");
      if (id.empty()) id = attr_value(t[8], "gene_id");
    }
    if (id.empty()) continue;
    // strip gene: prefix sometimes
    if (id.rfind("gene:", 0) == 0) id = id.substr(5);
    int64_t start = std::stoll(t[3]);
    int64_t end = std::stoll(t[4]);
    char strand = t[6].empty() ? '+' : t[6][0];
    int64_t tss = (strand == '-') ? end : start;
    if (out.count(id) || multi.count(id)) {
      if (out.count(id)) {
        warn("GFF multi-hit gene skipped: " + id);
        out.erase(id);
        multi.insert(id);
      }
      continue;
    }
    GeneLoc g;
    g.id = id;
    g.chrom = t[0];
    g.tss = tss;
    g.strand = strand;
    g.ok = true;
    out.emplace(id, std::move(g));
  }
  info("gff: " + std::to_string(n_gene) + " gene features, " +
       std::to_string(out.size()) + " unique ids with TSS");
  return out;
}

} // namespace eqtl

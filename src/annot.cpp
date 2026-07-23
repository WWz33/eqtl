#include "eqtl/annot.hpp"
#include "eqtl/util.hpp"

#include "gff3.hpp"

#include <unordered_set>

namespace eqtl {

static std::string first_attr(const std::string& attrs, const std::string& key) {
  auto parsed = gffsub::parse_attributes(attrs);
  auto it = parsed.find(key);
  if (it == parsed.end() || it->second.empty()) return {};
  return it->second.front();
}

static std::string gene_key(const gffsub::GffRecord& rec, const std::string& id_key) {
  if (!id_key.empty()) return first_attr(rec.attr_raw, id_key);
  if (rec.id && !rec.id->empty()) return *rec.id;
  std::string name = first_attr(rec.attr_raw, "Name");
  if (!name.empty()) return name;
  if (rec.gene_id && !rec.gene_id->empty()) return *rec.gene_id;
  return first_attr(rec.attr_raw, "gene_id");
}

std::unordered_map<std::string, GeneLoc> load_gff_tss(
    const std::string& path, const std::string& id_key) {
  gffsub::GffData data;
  gffsub::IdIndex idx;
  if (gffsub::parse_file(path, data, idx, gffsub::InputFormat::GFF3) != 0) {
    die("cannot parse gff via gffsub: " + path);
  }

  std::unordered_map<std::string, GeneLoc> out;
  std::unordered_set<std::string> multi;
  int n_gene = 0;

  for (const auto& rec : data) {
    if (rec.type != "gene") continue;
    ++n_gene;
    std::string id = gene_key(rec, id_key);
    if (id.empty()) continue;
    if (id.rfind("gene:", 0) == 0) id = id.substr(5);

    char strand = (rec.strand == '-' || rec.strand == '+') ? rec.strand : '+';
    int64_t tss = (strand == '-') ? rec.end : rec.start;

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
    g.chrom = rec.seqid;
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

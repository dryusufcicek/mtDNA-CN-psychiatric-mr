#!/usr/bin/env python
"""
11_instrument_vep_annotation.py
===============================
Functional annotation of the mtDNA-CN genetic instruments (Supplementary Table 18).
For every unique instrument rsID across the four exposure panels, query the Ensembl
Variant Effect Predictor REST endpoint (GRCh38.p14) and record the gene symbol, the
most-severe consequence and the PolyPhen / SIFT predictions. The Ensembl BioMart
"Human Short Variants" dataset returns identical PolyPhen/SIFT calls; the REST VEP
endpoint is used here for convenience.

Damaging = PolyPhen 'probably_damaging' OR SIFT 'deleterious'. This substantiates the
predicted-damaging missense instruments highlighted in the manuscript (DGUOK
rs74874677, POLG2 rs17850455, HABP2 rs7080536, PLAUR rs4760).

Reads:
  - data/instruments/unique_rsids.txt      (one rsID per line; union of the four panels)
Writes:
  - results/annotation/instrument_vep_annotation.tsv   (-> Supplementary Table 18)
"""
import os, json, time, urllib.request

BASE = os.environ.get("ITER003_BASE", ".")
IDS  = os.path.join(BASE, "data", "instruments", "unique_rsids.txt")
OUT  = os.path.join(BASE, "results", "annotation", "instrument_vep_annotation.tsv")
ENDPOINT = "https://rest.ensembl.org/vep/human/id"

# most-severe-consequence ranking (high -> low impact) for gene attribution
SEVERITY = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
            "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
            "missense_variant", "inframe_insertion", "inframe_deletion",
            "splice_region_variant", "synonymous_variant", "5_prime_UTR_variant",
            "3_prime_UTR_variant", "intron_variant", "upstream_gene_variant",
            "downstream_gene_variant", "regulatory_region_variant", "intergenic_variant"]


def rank(term):
    return SEVERITY.index(term) if term in SEVERITY else len(SEVERITY)


def vep_batch(ids):
    req = urllib.request.Request(
        ENDPOINT, data=json.dumps({"ids": ids}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def main():
    ids = [ln.strip() for ln in open(IDS) if ln.strip().startswith("rs")]
    out = {}
    for i in range(0, len(ids), 100):                      # VEP accepts up to 200 ids/POST
        for v in vep_batch(ids[i:i + 100]):
            rid = v.get("input"); msc = v.get("most_severe_consequence", "")
            best = None
            for tc in v.get("transcript_consequences", []):
                gs = tc.get("gene_symbol")
                if not gs:
                    continue
                tr = min([rank(t) for t in tc.get("consequence_terms", [])] + [len(SEVERITY)])
                score = (tr, 0 if tc.get("biotype") == "protein_coding" else 1)
                if best is None or score < best[0]:
                    best = (score, gs, tc.get("polyphen_prediction", ""), tc.get("sift_prediction", ""))
            gene, pph, sift = (best[1], best[2], best[3]) if best else ("", "", "")
            out[rid] = (gene, msc, pph, sift)
        time.sleep(1)                                      # be polite to the REST API

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write("rsID\tgene\tconsequence\tPolyPhen\tSIFT\n")
        for rid in ids:
            g, c, p, s = out.get(rid, ("", "not_found", "", ""))
            f.write(f"{rid}\t{g}\t{c}\t{p}\t{s}\n")
    dmg = sum(1 for g, c, p, s in out.values() if p == "probably_damaging" or s == "deleterious")
    print(f"annotated {len(out)}/{len(ids)} instruments; {dmg} predicted-damaging -> {OUT}")


if __name__ == "__main__":
    main()

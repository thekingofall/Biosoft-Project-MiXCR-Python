"""RepSeq.IO v1.6 reference exported from the verified MiXCR 3.0.13 artifact."""
from pathlib import Path
import json


class Reference:
    def __init__(self, document):
        self.metadata = {k: v for k, v in document.items() if k != "genes"}
        self.genes = document["genes"]
        self.used_gene_count = len(self.genes)
        self.selected_chains = "ALL"
        self.gene_hashes = {}
        self.by_name = {g["name"]: g for g in self.genes}
        self.by_type = {kind: [g for g in self.genes if g["type"] == kind]
                        for kind in "VDJC"}
        self._cache = {}
        for gene in self.genes:
            gene["chain_set"] = frozenset(gene["chains"].split(","))
            seq = gene["sequence"]
            gene["kmers"] = frozenset(seq[i:i+5] for i in range(len(seq)-4))

    def select_chains(self, chains):
        """Record the original aligner's reference selection for clone ties."""
        from .alignment import chain_names
        selected = chain_names(chains)
        self.used_gene_count = sum(bool(gene['chain_set'] & selected) for gene in self.genes)
        self.selected_chains = chains
        return self


def load_reference(path=None, v_feature='VRegionWithP'):
    if v_feature not in ('VRegionWithP','VGeneWithP'):
        raise ValueError('Supported V alignment features: VRegionWithP, VGeneWithP')
    default_name='human_reference.json' if v_feature=='VRegionWithP' else 'human_vgene_reference.json'
    path = Path(path) if path else Path(__file__).parent / "data" / default_name
    with path.open(encoding="utf-8") as stream:
        reference = Reference(json.load(stream))
    hash_path = path.with_name("human_reference_hashes.json")
    if hash_path.exists():
        reference.gene_hashes = json.loads(hash_path.read_text(encoding="utf-8"))
    else:
        reference.gene_hashes = {}
    feature_path = path.with_name("vgene_feature_points.json" if reference.metadata.get('v_feature')=='VGeneWithP' else "feature_points.json")
    reference.features = {}
    reference.point_descriptors = {}
    if feature_path.exists():
        metadata = json.loads(feature_path.read_text(encoding="utf-8"))
        reference.features = metadata["features"]
        reference.point_descriptors = metadata.get("points", {})
        for gene in reference.genes:
            gene["feature_points"] = metadata["genes"].get(gene["name"], {})
            gene["feature_intersections"] = metadata.get("intersections", {}).get(gene["name"], {})
    return reference

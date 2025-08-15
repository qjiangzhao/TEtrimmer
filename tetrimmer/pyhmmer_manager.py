# tetrimmer/pyhmmer_manager.py
from __future__ import annotations

import os
import io
import re
import sys
import time
import pickle
import collections
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterable, Any

import pyhmmer


# -------------------------------
# PFAM FILTERING / PRETTY PRINT
# -------------------------------

class pfam_filterer:
    def __init__(self, pfam_data_file: str | os.PathLike):
        self.pfam_file = str(pfam_data_file)
        self.pfam_data: Dict[str, Any] = {}
        # TE name like "rnd_1_family_470[261-1628]" or with trailing "_1"
        self.te_name_regex = re.compile(r'(.+?)\[(\-?\d+)\-(\-?\d+)](?:_\d+)?$')
        self.load_pfam_data()

    def load_pfam_data(self) -> None:
        """Reads the Pfam data file (Pfam-A.hmm.dat) to a dictionary."""
        Data = collections.namedtuple('Data', ['type', 'clan', 'ga_seq', 'ga_dom'])
        self.pfam_data = {}
        hmm_name = None
        typ = None
        clan = None
        ga_seq = ga_dom = None

        pfam_path = Path(self.pfam_file)
        if not pfam_path.is_file():
            raise FileNotFoundError(f"Pfam .dat file not found: {pfam_path}")

        with pfam_path.open("r") as fh:
            for line in fh:
                if line.startswith('#=GF ID'):
                    hmm_name = line[10:-1]
                elif line.startswith('#=GF TP'):
                    typ = line[10:-1]
                elif line.startswith('#=GF CL'):
                    clan = line[10:-1]
                elif line.startswith('#=GF GA'):
                    scores = line[10:-1].strip().rstrip(';').split(';')
                    ga_seq = float(scores[0])
                    ga_dom = float(scores[1])
                elif line.startswith('//'):
                    if hmm_name is not None:
                        self.pfam_data[hmm_name] = Data(typ, clan, ga_seq, ga_dom)
                    hmm_name = typ = clan = None
                    ga_seq = ga_dom = None

    def parse_hmmsearch_output(self, hmmsearch_records: Iterable[pyhmmer.plan7.TopHits]) -> Dict[str, List[Any]]:
        DOMAIN_INFO = [
            'te_name', 'te_start', 'te_end',
            'aln_start', 'aln_end',
            'domain_start', 'domain_end',
            'hmm_acc', 'hmm_name',
            'direction', 'evalue', 'clan'
        ]
        Domain = collections.namedtuple('Domain', DOMAIN_INFO)

        results: Dict[str, List[Any]] = {}
        for tophits in hmmsearch_records:
            hmm_accession = (tophits.query.accession or b"").decode("ascii", errors="ignore")
            hmm_name = (tophits.query.name or b"").decode("ascii", errors="ignore")
            # Some models might be missing in dat; guard it
            clan = self.pfam_data.get(hmm_name).clan if hmm_name in self.pfam_data else None

            for hit in tophits:
                name_long = (hit.name or b"").decode("ascii", errors="ignore")
                m = self.te_name_regex.search(name_long)
                if not m:
                    # Skip weirdly named sequences rather than crashing
                    continue
                seq_name, te_start_s, te_end_s = m.groups()
                te_start = int(te_start_s)
                te_end = int(te_end_s)
                aln = hit.best_domain.alignment
                align_from = int(aln.target_from)
                align_to = int(aln.target_to)
                evalue = float(hit.evalue)

                if te_end >= te_start:
                    direction = "+"
                    domain_start = te_start + (3 * align_from) - 3
                    domain_end = te_start + (3 * align_to) - 1
                else:
                    direction = "-"
                    domain_start = te_start - (3 * align_to) + 1
                    domain_end = te_start - (3 * align_from) + 3

                dom = Domain(
                    seq_name, te_start, te_end,
                    align_from, align_to,
                    domain_start, domain_end,
                    hmm_accession, hmm_name,
                    direction, evalue, clan
                )
                results.setdefault(seq_name, []).append(dom)

        # Sort domains by alignment start per sequence
        for seq_name in results:
            results[seq_name].sort(key=lambda x: x.aln_start)

        return results

    def resolve_overlapping_domains(self, results: Dict[str, List[Any]]) -> Dict[str, List[Any]]:
        """Resolves overlapping domains belonging to the same clan."""
        for seq_id, domains in list(results.items()):
            changed = True
            while changed and len(domains) > 1:
                changed = False
                to_drop = set()
                for i in range(len(domains) - 1):
                    di = domains[i]
                    dj = domains[i + 1]
                    if di.clan and dj.clan and di.clan == dj.clan and dj.aln_start <= di.aln_end:
                        # Keep the better (lower e-value)
                        drop_idx = i if di.evalue > dj.evalue else i + 1
                        to_drop.add(drop_idx)
                        changed = True
                if changed:
                    domains = [d for idx, d in enumerate(domains) if idx not in to_drop]
            results[seq_id] = domains
        return results

    def pretty_print(self, results: Dict[str, List[Any]], outfile: str | os.PathLike) -> str:
        """Write results as TSV. Ensures parent dir exists. Returns path string."""
        if not outfile:
            raise ValueError("pretty_print: outfile must be a non-empty path string")

        of = Path(outfile)
        of.parent.mkdir(parents=True, exist_ok=True)

        header = ['TE_name', 'orf_start', 'orf_end', 'domain_start', 'domain_end', 'direction', 'domain_name', 'domain_reference']
        with of.open("w") as out:
            print(*header, sep="\t", file=out)
            for seq in sorted(results):
                for domain in results[seq]:
                    if domain.direction == "+":
                        print(seq, domain.te_start, domain.te_end,
                              domain.domain_start, domain.domain_end,
                              domain.direction, domain.hmm_name, domain.hmm_acc,
                              sep="\t", file=out)
                    else:
                        print(seq, domain.te_end, domain.te_start,
                              domain.domain_start, domain.domain_end,
                              domain.direction, domain.hmm_name, domain.hmm_acc,
                              sep="\t", file=out)
        return str(of)


# -------------------------------
# HMMER MANAGER
# -------------------------------

class pyhmmer_manager:
    def __init__(self, hmm_models: str | os.PathLike | None = None,
                 hmm_dat: str | os.PathLike | None = None,
                 cpus: int = 1,
                 default_outdir: str | os.PathLike | None = None):
        self.model_file = str(hmm_models) if hmm_models is not None else None
        self.dat_file = str(hmm_dat) if hmm_dat is not None else None
        self.threads = max(1, int(cpus))
        self.default_outdir = Path(default_outdir) if default_outdir else Path.cwd() / "pfam_out"

        self.hmm_models: Optional[List[pyhmmer.plan7.OptimizedProfile]] = None
        self.nt_alphabet = pyhmmer.easel.Alphabet.dna()
        self.aa_alphabet = pyhmmer.easel.Alphabet.amino()
        self.pfam_operator: Optional[pfam_filterer] = None

    # ---- Preparation

    def prep_pfam_op(self) -> None:
        if not self.dat_file:
            raise ValueError("Pfam .dat path is not set")
        self.pfam_operator = pfam_filterer(self.dat_file)

    def load_hmm_database(self) -> None:
        """Load HMMs and convert to optimized profiles."""
        if not self.model_file:
            raise ValueError("HMM models path is not set")

        hmm_path = Path(self.model_file)
        if not hmm_path.is_file():
            raise FileNotFoundError(f"HMM file not found: {hmm_path}")

        print("\nLoading DB")
        print("Interpreting models")
        # Use path directly; no need to slurp entire file in memory
        with pyhmmer.plan7.HMMFile(str(hmm_path)) as fh:
            hmms = list(fh)

        print("Optimizing HMMs for search")
        # to_profile().to_optimized() is the fast path for searches
        self.hmm_models = [h.to_profile().to_optimized() for h in hmms]

        print("HMM search prepared!\n")

    def prepare(self) -> None:
        self.prep_pfam_op()
        self.load_hmm_database()

    # ---- Sequence preparation

    def _read_fasta_as_aa(self, fasta_path: str | os.PathLike) -> Tuple[List[pyhmmer.easel.DigitalSequence], int]:
        """Reads a FASTA (expected to be amino-acid ORFs) and returns digitized sequences."""
        p = Path(fasta_path)
        if not p.is_file():
            raise FileNotFoundError(f"Input FASTA not found: {p}")

        total_len = 0
        seqs: List[pyhmmer.easel.DigitalSequence] = []
        # Use Easel's FASTA reader for robustness
        with pyhmmer.easel.SequenceFile(str(p), digital=True, alphabet=self.aa_alphabet) as sf:
            for dsq in sf:
                total_len += len(dsq)
                # Ensure sequence has a name
                if not dsq.name:
                    dsq.name = p.name.encode("ascii", errors="ignore")
                seqs.append(dsq)
        return seqs, total_len

    # Backwards-compatible name (your original API)
    def prepare_nucleotide_seq_for_search(self, nucleotide_fasta: str | os.PathLike):
        # Many TEtrimmer inputs here are actually ORF AA sequences; keep AA default.
        return self._read_fasta_as_aa(nucleotide_fasta)

    # ---- Search

    def execute_search(self,
                       sequences: Iterable[pyhmmer.easel.DigitalSequence],
                       seq_e_value_cutoff: float = 1e-2,
                       domain_e_value_cutoff: float = 1e-2,
                       output_file: Optional[str | os.PathLike] = None,
                       cutoff_method: str = "evalue") -> bool:
        """Run hmmsearch and pretty-print PFAM domains to output_file."""
        if not self.hmm_models:
            raise RuntimeError("HMM models not loaded. Call prepare() first.")
        if not self.pfam_operator:
            raise RuntimeError("PFAM operator not prepared. Call prepare() first.")

        # Choose cutoff mode
        cutoff_method = cutoff_method.lower().strip()
        if cutoff_method == "evalue":
            res = pyhmmer.hmmsearch(
                self.hmm_models, sequences,
                cpus=self.threads, E=seq_e_value_cutoff, domE=domain_e_value_cutoff,
                Z=len(self.hmm_models)
            )
        elif cutoff_method == "gathering":
            res = pyhmmer.hmmsearch(
                self.hmm_models, sequences,
                bit_cutoffs="gathering", cpus=self.threads, Z=len(self.hmm_models)
            )
        elif cutoff_method == "trusted":
            res = pyhmmer.hmmsearch(
                self.hmm_models, sequences,
                bit_cutoffs="trusted", cpus=self.threads, Z=len(self.hmm_models)
            )
        else:
            raise ValueError(f"Unknown cutoff_method: {cutoff_method}")

        # Keep only models with at least one hit
        res = [r for r in res if len(r) > 0]

        formatted = self.pfam_operator.parse_hmmsearch_output(res)
        cleaned = self.pfam_operator.resolve_overlapping_domains(formatted)

        # Determine if any domains were found
        has_any_domains = any(len(domains) > 0 for domains in cleaned.values())

        # Guarantee an output path
        outfile = self._ensure_outfile(output_file)
        self.pfam_operator.pretty_print(cleaned, outfile)

        return has_any_domains

    def _ensure_outfile(self, output_file: Optional[str | os.PathLike]) -> str:
        """Ensure we have a valid path to write; auto-generate if None."""
        if output_file:
            of = Path(output_file)
            of.parent.mkdir(parents=True, exist_ok=True)
            return str(of)
        # Fallback path
        self.default_outdir.mkdir(parents=True, exist_ok=True)
        fname = f"pfam_{os.getpid()}_{int(time.time()*1000)}.tsv"
        return str(self.default_outdir / fname)

    # ---- Runner (used by your process handler)

    def run(self, args) -> Tuple[str, bool]:
        """
        Accepts either:
          - dict: {'infile': <path>, 'outfile': <path or None>}
          - tuple/list: (infile, outfile)
        Returns: (outfile_path, domains_detected)
        """
        if isinstance(args, dict):
            infile = args.get("infile") or args.get("seq_file") or args.get("input") or args.get("path")
            outfile = args.get("outfile")
        else:
            try:
                infile, outfile = args[0], args[1]
            except Exception:
                raise TypeError("pyhmmer_manager.run expects dict or (infile, outfile)")

        if infile is None or isinstance(infile, bool):
            raise TypeError(f"infile must be a path string, got: {type(infile).__name__}={infile}")
        if outfile is not None and isinstance(outfile, bool):
            # Prevent the "bool where path expected" error later
            raise TypeError(f"outfile must be path-like or None, got bool: {outfile}")

        seqs, _ = self.prepare_nucleotide_seq_for_search(infile)
        outfile_path = self._ensure_outfile(outfile)
        domains_detected = self.execute_search(sequences=seqs, output_file=outfile_path)
        return (outfile_path, domains_detected)

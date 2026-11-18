#!/usr/bin/env python3
"""Reconstruct a graph-supported circular Flye candidate without inventing bases.

NANOPlasty is used to select an exactly matching circular traversal; therefore
this candidate is not independent confirmation of plasmid circularity.
"""
from itertools import permutations, product
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def rc(seq):
    return str(Seq(seq).reverse_complement())


def equivalence(a, b):
    return len(a) == len(b) and (a in b + b or a in rc(b) + rc(b))


def longest_overlap(a, b, minimum=20, maximum=500):
    for n in range(min(len(a), len(b), maximum), minimum - 1, -1):
        if a[-n:] == b[:n]:
            return n
    return 0


def graph_data(path):
    links, paths = set(), {}
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] == "L" and len(f) >= 6:
                a, sa, b, sb = f[1:5]
                links.add((a, sa, b, sb))
                flip = {"+": "-", "-": "+"}
                links.add((b, flip[sb], a, flip[sa]))
            elif f[0] == "P" and len(f) >= 3:
                steps = f[2].split(",")
                if len(steps) == 1 and steps[0][-1] in "+-":
                    paths[f[1]] = (steps[0][:-1], steps[0][-1])
    return links, paths


def reconstruct(contigs, graph_path, nano_path, expected_length):
    if not 2 <= len(contigs) <= 6:
        raise ValueError("Graph reconstruction supports 2–6 contigs")
    with open(nano_path) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) != 1:
        raise ValueError("Expected one NANOPlasty contig")
    nano = str(records[0].seq).upper()
    links, paths = graph_data(graph_path)
    names = sorted(contigs)
    if any(name not in paths for name in names):
        raise ValueError("Flye graph lacks a one-edge path for a contig")
    candidates = []
    # Fix one contig to avoid enumerating rotations of the same cycle.
    for tail in permutations(names[1:]):
        order = (names[0],) + tail
        for strands in product("+-", repeat=len(order)):
            nodes = []
            for name, direction in zip(order, strands):
                edge, strand = paths[name]
                strand = strand if direction == "+" else ("-" if strand == "+" else "+")
                nodes.append((edge, strand))
            if any((nodes[i][0], nodes[i][1],
                    nodes[(i+1) % len(nodes)][0], nodes[(i+1) % len(nodes)][1])
                   not in links for i in range(len(nodes))):
                continue
            sequences = [contigs[name] if strand == "+" else rc(contigs[name])
                         for name, strand in zip(order, strands)]
            overlaps = [longest_overlap(sequences[i], sequences[(i+1) % len(order)])
                        for i in range(len(order))]
            # Trim the duplicate prefix of the downstream segment at each junction.
            candidate = "".join(sequences[i][overlaps[i-1]:]
                                for i in range(len(order)))
            if not 0.7 * expected_length <= len(candidate) <= 1.3 * expected_length:
                continue
            if equivalence(candidate, nano):
                candidates.append((candidate, order, strands, overlaps))
    if not candidates:
        raise ValueError("No graph-connected traversal matches NANOPlasty circularly; "
                         "raw Flye contigs preserved for manual review")
    first = candidates[0]
    if any(not equivalence(x[0], first[0]) for x in candidates):
        raise ValueError("Several non-equivalent candidate plasmid sequences")
    return first, len(candidates)


def process(args, write_fasta, write_tsv):
    with open(args.assembly) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) < 2:
        raise ValueError("Graph reconstruction requires multiple contigs")
    contigs = {r.id: str(r.seq).upper() for r in records}
    if len(contigs) != len(records):
        raise ValueError("Repeated contig names")
    if any(not seq or set(seq) - set("ACGTN") for seq in contigs.values()):
        raise ValueError("Unexpected Flye sequence bases")
    with open(args.design) as handle:
        design = list(SeqIO.parse(handle, "fasta"))
    if len(design) != 1:
        raise ValueError("Expected one designed plasmid FASTA")
    (seq, order, strands, overlaps), count = reconstruct(
        contigs, args.graph, args.nanoplasty, len(design[0].seq))
    sample = args.sample or "flye"
    for path in (args.output_fasta, args.report, args.boundaries,
                 args.pairwise, args.copies_fasta):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.output_fasta, [(sample + "_flye_graph_reconstructed", seq)])
    write_fasta(args.copies_fasta, [(name, contigs[name]) for name in sorted(contigs)])
    write_tsv(args.boundaries, [
        {"copy": i+1, "start_0based": "", "end_0based_exclusive": "",
         "wraps_raw_origin": "not_applicable", "length_bp": len(contigs[name])}
        for i, name in enumerate(order)],
        ["copy", "start_0based", "end_0based_exclusive",
         "wraps_raw_origin", "length_bp"])
    write_tsv(args.pairwise, [],
              ["copy_a", "copy_b", "identities", "mismatches", "gap_bases",
               "identity", "differences", "aligned_columns"])
    write_tsv(args.report, [{
        "sample": sample,
        "status": "graph_reconstructed_nanoplasty_concordant_read_review_required",
        "raw_contig_id": ",".join(order),
        "raw_length_bp": sum(len(s) for s in contigs.values()),
        "design_length_bp": len(design[0].seq),
        "estimated_copies": "not_applicable",
        "anchor_length_bp": "not_applicable",
        "anchor_reference_start_0based": "not_applicable",
        "raw_orientation_for_detection": ",".join(strands),
        "copy_lengths_bp": ",".join(str(len(contigs[n])) for n in order),
        "selected_copy": "not_applicable",
        "normalized_length_bp": len(seq),
        "min_pairwise_identity": "not_applicable",
        "notes": "graph path " + ",".join(n+s for n,s in zip(order,strands))
                 + "; overlaps-to-next " + ",".join(map(str,overlaps))
                 + f"; matching traversals {count}; exact circular equivalence to NANOPlasty; "
                   "raw-read junction review still required"
    }], ["sample", "status", "raw_contig_id", "raw_length_bp", "design_length_bp",
         "estimated_copies", "anchor_length_bp", "anchor_reference_start_0based",
         "raw_orientation_for_detection", "copy_lengths_bp", "selected_copy",
         "normalized_length_bp", "min_pairwise_identity", "notes"])
    print(f"{sample}: graph-reconstructed {len(seq)} bp; "
          f"path={','.join(n+s for n,s in zip(order,strands))}; "
          f"overlaps={overlaps}; read review required", flush=True)

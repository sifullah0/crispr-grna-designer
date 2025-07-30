import re
import RNA

def find_forward_grna_candidates(target_sequence, pam_pattern=r'[ACGT]GG'):
    """
    Scans the positive‑sense DNA for candidate gRNAs:
    - 23‑nt windows ending in NGG
    - Extracts the 20‑nt spacer upstream of the PAM
    Returns a list of dicts: {"sequence", "start", "end"}
    """
    seq = target_sequence.upper()
    candidates = []

    for i in range(len(seq) - 22):
        window = seq[i : i + 23]
        spacer, pam = window[:20], window[20:]
        if re.fullmatch(pam_pattern, pam):
            candidates.append({
                "sequence": spacer,
                "start": i + 1,          # 1‑based
                "end": i + 23            # inclusive end
            })
    return candidates


def gRNA_score(gRNA_seq):
    """
    Takes a 20‑nt RNA guide, folds it, and returns either:
     - ('mfe', mfe_value, structure)   if mfe < -15
     - (avg_unpaired_prob, structure)   otherwise
    """
    fc = RNA.fold_compound(gRNA_seq)
    structure, mfe = fc.mfe()

    if mfe < -15:
        return ('mfe', mfe, structure)
    else:
        fc.exp_params_rescale(mfe)
        fc.pf()
        plist = fc.plist_from_probs(1e-6)
        bpp = {(ep.i, ep.j): ep.p for ep in plist}

        unpaired_probs = []
        for i in range(1, min(13, len(gRNA_seq) + 1)):
            paired = sum(p for (a, b), p in bpp.items() if a == i or b == i)
            unpaired_probs.append(1.0 - paired)

        avg_unpaired = sum(unpaired_probs) / len(unpaired_probs)
        return (avg_unpaired, structure)


if __name__ == "__main__":
    # Example target DNA
    long_dna = (
        "ATGGGCCTACGTTAGCTAGCTAGCTAGGCTAATTCGGAATCGATCGATGCTAGCTGATCG"
        "GTTAGCTAGGTACGATCGTAGCTAGCTAGCTAGCTAGGCTAACGATCGTAGCTAGCTAAG"
        "CTAGGCTAGCTAGCTAACGTTAGCTAGCAGCTAGCTAGCTAGCTAACGATCGTAGCTAGC"
    )

    # 1. Extract forward‑strand candidates
    candidates = find_forward_grna_candidates(long_dna)

    # 2. Score each candidate
    scored_list = []
    for cand in candidates:
        seq = cand["sequence"]
        result = gRNA_score(seq)
        # unify the return format: (type/score, structure, optional MFE)
        if result[0] == 'mfe':
            _, mfe_value, struct = result
            score = mfe_value
        else:
            score, struct = result
        scored_list.append({
            "sequence": seq,
            "start": cand["start"],
            "end": cand["end"],
            "score": score,
            "structure": struct
        })

    # 3. Sort by score descending (higher is better)
    scored_list.sort(key=lambda x: x["score"], reverse=True)

    # 4. Print top candidates
    print(f"Found {len(scored_list)} candidates, sorted by score:\n")
    for entry in scored_list:
        print(f"{entry['sequence']} ({entry['start']}-{entry['end']}) | "
              f"Score: {entry['score']:.3f} | Structure: {entry['structure']}")

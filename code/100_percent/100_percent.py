import re
import RNA
import heapq
import azimuth.model_comparison
import pandas as pd
import numpy as np
import warnings

# Suppress all warnings
warnings.filterwarnings("ignore")

def find_forward_grna_candidates(target_sequence, pam_pattern=r'[ACGT]GG'):
    """
    Scans the positive-sense DNA for candidate gRNAs:
    - 23-nt windows ending in NGG
    - Extracts the 20-nt spacer upstream of the PAM
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
                "start": i + 1,          # 1-based
                "end": i + 23            # inclusive end
            })
    return candidates

def filter_grna(gRNA_seq):
    """
    Filters gRNAs based on:
    - GC content between 40% and 60%
    - Absence of poly-T motif ("TTTT")
    Returns True if gRNA passes filters, False otherwise
    """
    gc_count = gRNA_seq.count('G') + gRNA_seq.count('C')
    gc_content = gc_count / len(gRNA_seq)
    if not (0.4 <= gc_content <= 0.6):
        return False
    if 'TTTT' in gRNA_seq:
        return False
    return True

def gRNA_score(gRNA_seq):
    """
    Takes a 20-nt RNA guide, folds it, and returns either:
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
        return (avg_unpaired, structure, mfe)

def get_azimuth_score(gRNA_seq, target_sequence, start):
    seq = target_sequence.upper()
    if start < 5 or start > len(seq) - 25:
        return 0.0
    sequence_30mer = seq[start - 5:start + 25]  # 1-based adjustment
    if len(sequence_30mer) != 30:
        return 0.0
    data = pd.DataFrame({'sequence': [sequence_30mer]})
    scores = azimuth.model_comparison.predict(data['sequence'].to_numpy())
    return scores[0]

if __name__ == "__main__":
    # Example target DNA
    long_dna = (
        "ATGGGCCTACGTTAGCTAGCTAGCTAGGCTAATTCGGAATCGATCGATGCTAGCTGATCG"
        "GTTAGCTAGGTACGATCGTAGCTAGCTAGCTAGCTAGGCTAACGATCGTAGCTAGCTAAG"
        "CTAGGCTAGCTAGCTAACGTTAGCTAGCAGCTAGCTAGCTAGCTAACGATCGTAGCTAGC"
        "AGTCTGATCGAGTCCGAGCAGAAGAGGATAGAGTCTGATCGAGTCCGAGCAGAAGGGATA"
    )

    # 1. Extract forward-strand candidates
    candidates = find_forward_grna_candidates(long_dna)

    # 2. Filter candidates by GC content and poly-T
    filtered_candidates = [cand for cand in candidates if filter_grna(cand["sequence"])]

    # 3. Score each candidate
    scored_list = []
    heap = []
    for cand in filtered_candidates:
        seq = cand["sequence"]
        # Get Azimuth score
        azimuth_score = get_azimuth_score(seq, long_dna, cand["start"])
        # Get accessibility score
        result = gRNA_score(seq)
        if result[0] == 'mfe':
            _, mfe_value, struct = result
            accessibility_score = mfe_value
        else:
            accessibility_score, struct, mfe = result
        # Calculate total score: 0.6 * Azimuth + 0.4 * Accessibility
        total_score = 0.6 * azimuth_score + 0.4 * accessibility_score
        # Push to heap for sorting (use negative score for max-heap)
        heapq.heappush(heap, (-total_score, seq, cand["start"], cand["end"], struct, accessibility_score, azimuth_score, mfe))

    # 4. Extract sorted candidates from heap
    while heap:
        neg_score, seq, start, end, struct, accessibility_score, azimuth_score, mfe = heapq.heappop(heap)
        scored_list.append({
            "sequence": seq,
            "start": start,
            "end": end,
            "accessibility_score": accessibility_score,
            "azimuth_score": azimuth_score,
            "total_score": -neg_score,  # Store total score
            "mfe": mfe,
            "structure": struct
        })

    # 5. Print top candidates
    print(f"Found {len(scored_list)} candidates, sorted by total score:\n")
    for entry in scored_list:
        print(f"{entry['sequence']} ({entry['start']}-{entry['end']}) | "
              f"Azimuth: {entry['azimuth_score']:.3f} | "
              f"Accessibility: {entry['accessibility_score']:.3f} | "
              f"Total Score: {entry['total_score']:.3f} | MFE: {entry['mfe']:.2f} | Structure: {entry['structure']}")
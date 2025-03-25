import numpy as np
from math import log

def invasion_score_from_counts(motif_data, neutral_score=1.0, max_sites=500):
    """
    Adds confidence weighting based on motif site counts.
    """
    scores = []
    weights = []
    total_sites = 0

    for m in motif_data:
        h_total = m['host_total']
        h_meth = m['host_meth']
        p_total = m['plasmid_total']
        p_meth = m['plasmid_meth']

        if h_total == 0:
            continue

        weight = h_total
        total_sites += h_total + p_total

        f_host = h_meth / h_total
        if p_total == 0:
            motif_score = neutral_score
        else:
            f_plasmid = p_meth / p_total
            motif_score = 1 - abs(f_host - f_plasmid)

        scores.append(motif_score * weight)
        weights.append(weight)

    if not scores:
        return {'invasion_score': 0.0, 'confidence': 0.0, 'final_score': 0.0}

    invasion_score = sum(scores) / sum(weights)

    # Confidence scaling (logarithmic)
    confidence = log(1 + total_sites) / log(1 + max_sites)
    final_score = invasion_score * confidence

    return {
        'invasion_score': round(invasion_score, 4),
        'confidence': round(confidence, 4),
        'final_score': round(final_score, 4),
        'total_sites': total_sites
    }




motif_data = [
    {'motif': 'GATC', 'host_total': 100, 'host_meth': 90, 'plasmid_total': 80, 'plasmid_meth': 70},
    {'motif': 'TTAA', 'host_total': 50, 'host_meth': 45, 'plasmid_total': 0, 'plasmid_meth': 0},
    {'motif': 'CCWGG', 'host_total': 20, 'host_meth': 18, 'plasmid_total': 20, 'plasmid_meth': 5},
]

result = invasion_score_from_counts(motif_data)
print(result)

# Output:
# {
#   'invasion_score': 0.8738,
#   'confidence': 0.6845,
#   'final_score': 0.5985,
#   'total_sites': 270
# }

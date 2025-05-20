import numpy as np
from scipy.special import betaln, logsumexp

def compute_p_same_room(x_data, y_data, alpha=1, beta=1):
    """
    Computes the posterior probability that x and y went to the same room.

    Parameters:
    - x_data: list of tuples (k_i, n_i) for person x, one per box
    - y_data: list of tuples (k_i, n_i) for person y, one per box
    - alpha, beta: prior parameters for Beta distribution (default uniform)

    Returns:
    - p_same: posterior probability that x and y went to the same room
    """

    def approx_log_beta_binomial(k, n, alpha, beta):
        # Approximate log marginal likelihood (ignores binomial coefficient)
        return betaln(k + alpha, n - k + beta) - betaln(alpha, beta)

    # Likelihood under same room: combine counts
    log_A = 0
    for (kx, nx), (ky, ny) in zip(x_data, y_data):
        k = kx + ky
        n = nx + ny
        log_A += approx_log_beta_binomial(k, n, alpha, beta)

    # Likelihood under different rooms: compute separately and sum
    log_B = 0
    for (kx, nx) in x_data:
        log_B += approx_log_beta_binomial(kx, nx, alpha, beta)
    for (ky, ny) in y_data:
        log_B += approx_log_beta_binomial(ky, ny, alpha, beta)

    # Compute posterior
    log_total = logsumexp([log_A, log_B])
    log_p_same = log_A - log_total
    p_same = np.exp(log_p_same)

    return p_same


if __name__ == "__main__":
    # # Example with 2 boxes:
    # Dx = [(496, 496), (1439, 2329)]  # (red, total) for person x
    # Dy = [(0, 0), (3, 9)]            # (red, total) for person y

    Dx = [(33813,41250), (3136, 3405)]  # (red, total) for person x
    Dy = [(322, 324), (104, 104)]            # (red, total) for person y

    p_same = compute_p_same_room(Dx, Dy)
    print(f"P(same room) = {p_same:.4f}")

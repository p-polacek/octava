#!/usr/bin/env python
import os
import json
import numpy as np

from app.dirs import data_dir

def kmed(genome, markers, N, include=[], exclude=[], max_iter=100):
    kin_fn = os.path.join(data_dir, genome, "assoc", markers, "stats", "kin")

    with open(kin_fn) as kin_in:
        k, labels = json.load(kin_in)

    k_matr = np.array(k)
    d_matr = np.max(k_matr) - k_matr

    if any(samp not in labels for samp in include):
        return [[], "Some pre-set samples are missing from selected Marker set."]
    if set(include).intersection(set(exclude)):
        return [[], "Some samples are selected to both include and exclude."]
    if len(set(labels).difference(set(exclude))) <= N:
        return [[], "After removing samples to exclude, remaining sample set must be larger than requested core set size."]

    if exclude:
        keep_idx = [labels.index(samp) for samp in labels if samp not in exclude]
        labels = [samp for samp in labels if samp not in exclude]
        d_matr = d_matr[keep_idx, :][:, keep_idx]

    n_samp = d_matr.shape[0]
    np.fill_diagonal(d_matr, 0)

    include = [labels.index(samp) for samp in include]
    selected = include[:]

    remaining = list(set(range(n_samp)) - set(selected))
    nearest_cost = np.full(n_samp, np.inf)

    # greedy fill (distance-based)
    while len(selected) < N:
        best = None
        best_score = np.inf

        for i in remaining:
            # Only need to consider what *adding i* contributes on top of current assignment
            updated = np.minimum(nearest_cost, d_matr[:, i])
            score = updated.sum()
            if score < best_score:
                best_score = score
                best = i

        selected.append(best)
        remaining.remove(best)

        # Update nearest costs with the newly added medoid
        nearest_cost = np.minimum(nearest_cost, d_matr[:, best])

    selected = np.array(selected)

    def compute_cost(medoids):
        dist = d_matr[:, medoids]
        return np.sum(np.min(dist, axis=1))

    # Precompute nearest and second-nearest medoid distances for each point
    def get_nearest(selected):
        dist = d_matr[:, selected]
        sorted_idx = np.argsort(dist, axis=1)
        nearest = dist[np.arange(n_samp), sorted_idx[:, 0]]
        second  = dist[np.arange(n_samp), sorted_idx[:, 1]] if len(selected) > 1 else np.full(n_samp, np.inf)
        nearest_medoid_pos = sorted_idx[:, 0]  # position in selected[], not sample index
        return nearest, second, nearest_medoid_pos

    selected = list(selected)
    include = set(include)
    selected_set = set(selected)

    current_cost = compute_cost(selected)

    for q in range(max_iter):
        nearest, second_nearest, nearest_medoid_pos = get_nearest(selected)

        best_delta = 0
        best_swap  = None

        non_medoid = np.array([j for j in range(n_samp) if j not in selected_set])
        D_cands = d_matr[:, non_medoid]  # (n_samp, n_candidates) — compute once per iteration

        for i, m in enumerate(selected):
            if m in include:
                continue

            is_nearest_m = nearest_medoid_pos == i  # (n_samp,)

            # For each candidate o (column), compute per-point delta vectorised
            # Shape: (n_samp, n_candidates)
            deltas = np.where(
                ~is_nearest_m[:, None],
                np.minimum(nearest[:, None], D_cands) - nearest[:, None],
                np.minimum(second_nearest[:, None], D_cands) - nearest[:, None],
            )

            total_deltas = deltas.sum(axis=0)  # (n_candidates,)
            best_j_pos = np.argmin(total_deltas)

            if total_deltas[best_j_pos] < best_delta:
                best_delta = total_deltas[best_j_pos]
                best_swap  = (i, non_medoid[best_j_pos])

        if best_swap is None:
            break

        i, j = best_swap
        selected_set.discard(selected[i])
        selected[i] = j
        selected_set.add(j)
        current_cost += best_delta

    core_set = [labels[samp] for samp in sorted(selected)]
    return [core_set, None]

#kmed = KMedoids(n_clusters=N, metric="precomputed", init="k-medoids++")#, random_state=0)
#kmed.fit_predict(d)
#selected = kmed.medoid_indices_

#def kmedoids_cost(D, medoids):
#    medoids = np.array(medoids)
#    dist = D[:, medoids]
#    return np.sum(np.min(dist, axis=1))

#cost_custom = kmedoids_cost(d, custom)
#cost_sklearn = kmedoids_cost(d, selected)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bayesian surrogate (GP) for sequence-level candidate prioritization.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel as C
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel
except ImportError:  # pragma: no cover - optional dependency
    GaussianProcessRegressor = None
    C = RBF = WhiteKernel = None


class BOSurrogateGP:
    """Simple GP surrogate with one-hot sequence encoding."""

    AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
    AA_TO_IDX = {aa: i for i, aa in enumerate(AA_ALPHABET)}

    def __init__(self):
        if GaussianProcessRegressor is None:
            raise ImportError(
                "scikit-learn is required for BO surrogate. "
                "Install scikit-learn or disable bo_enabled."
            )

        kernel = (
            C(1.0, (1e-3, 1e3))
            * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2))
            + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-6, 1e1))
        )
        self.model = GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=2,
            random_state=42,
        )
        self._observations: Dict[str, float] = {}
        self._sequence_len = None
        self._is_fitted = False

    @property
    def n_observations(self) -> int:
        return len(self._observations)

    def add_observation(self, sequence: str, target: float):
        sequence = str(sequence).strip().upper()
        if not sequence:
            return

        if self._sequence_len is None:
            self._sequence_len = len(sequence)

        if len(sequence) != self._sequence_len:
            return

        self._observations[sequence] = float(target)

    def is_ready(self, min_train: int) -> bool:
        return (
            self._is_fitted
            and self._sequence_len is not None
            and self.n_observations >= int(min_train)
        )

    def fit_if_ready(self, min_train: int) -> bool:
        if self.n_observations < int(min_train):
            self._is_fitted = False
            return False

        X, y = self._build_training_matrix()
        self.model.fit(X, y)
        self._is_fitted = True
        return True

    def predict_mu_sigma(self, sequences: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
        if not self._is_fitted:
            raise RuntimeError("Surrogate model is not fitted.")

        X = self._encode_sequences(sequences)
        mu, sigma = self.model.predict(X, return_std=True)
        sigma = np.nan_to_num(sigma, nan=0.0, posinf=0.0, neginf=0.0)
        return mu, sigma

    def acquisition(self, sequences: Sequence[str], beta: float) -> np.ndarray:
        """LCB for minimization: lower is better."""
        mu, sigma = self.predict_mu_sigma(sequences)
        return mu - float(beta) * sigma

    def _build_training_matrix(self) -> Tuple[np.ndarray, np.ndarray]:
        sequences = list(self._observations.keys())
        X = self._encode_sequences(sequences)
        y = np.array([self._observations[s] for s in sequences], dtype=np.float64)
        return X, y

    def _encode_sequences(self, sequences: Sequence[str]) -> np.ndarray:
        rows: List[np.ndarray] = []
        for seq in sequences:
            seq = str(seq).strip().upper()
            if self._sequence_len is None:
                self._sequence_len = len(seq)
            if len(seq) != self._sequence_len:
                raise ValueError(
                    f"Inconsistent sequence length for surrogate: "
                    f"{len(seq)} != {self._sequence_len}"
                )

            row = np.zeros(self._sequence_len * len(self.AA_ALPHABET), dtype=np.float32)
            for i, aa in enumerate(seq):
                aa_idx = self.AA_TO_IDX.get(aa)
                if aa_idx is None:
                    continue
                row[i * len(self.AA_ALPHABET) + aa_idx] = 1.0
            rows.append(row)

        return np.vstack(rows) if rows else np.empty((0, 0), dtype=np.float32)

import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, EsmForMaskedLM


class ESM2ProbMatrix:
    def __init__(self, model_name="facebook/esm2_t6_8M_UR50D"):
        self.device = torch.device(
            "cuda" if torch.cuda.is_available()
            else ("mps" if torch.backends.mps.is_available() else "cpu")
        )

        self.tokenizer = AutoTokenizer.from_pretrained(
            model_name, do_lower_case=False
        )

        self.model = EsmForMaskedLM.from_pretrained(model_name).to(self.device)
        self.model.eval()

        self.AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
        self.aa_ids = torch.tensor(
            [self.tokenizer.convert_tokens_to_ids(a) for a in self.AMINO_ACIDS],
            device=self.device
        )

        # Cache: seq -> dict(ll, log_probs, prob_matrix)
        self.cache = {}

    @torch.no_grad()
    def _forward(self, seq):
        """Single forward pass, cached."""
        if seq in self.cache:
            return self.cache[seq]

        tokens = self.tokenizer(seq, return_tensors="pt").to(self.device)
        input_ids = tokens["input_ids"]

        outputs = self.model(**tokens)
        logits = outputs.logits[:, 1:-1, :]  # remove <cls> and <eos>

        log_probs = F.log_softmax(logits, dim=-1)
        probs = torch.exp(log_probs)

        # Extract probs only for the 20 canonical AAs
        prob_matrix = probs[0, :, self.aa_ids]      # [L, 20]
        log_prob_matrix = log_probs[0, :, self.aa_ids]

        # True sequence likelihood
        true_ids = input_ids[0, 1:-1]
        ll_positions = log_probs[0, torch.arange(len(seq)), true_ids]
        ll_total = ll_positions.sum().item()

        data = {
            "ll": ll_total,
            "ll_positions": ll_positions,
            "prob_matrix": prob_matrix
        }

        self.cache[seq] = data
        return data

    def get_esm_ll(self, seq):
        """Total log-likelihood."""
        return self._forward(seq)["ll"]

    def get_probability_matrix(self, seq):
        """L x 20 probability matrix."""
        return self._forward(seq)["prob_matrix"]

    def most_probable_replacement(self, seq, position):
        """
        Returns:
          (best_aa, probability)
        """
        data = self._forward(seq)
        probs = data["prob_matrix"][position].clone()

        wt_aa = seq[position]
        wt_idx = self.AMINO_ACIDS.index(wt_aa)
        probs[wt_idx] = 0.0

        best_idx = probs.argmax().item()
        return self.AMINO_ACIDS[best_idx], probs[best_idx].item()

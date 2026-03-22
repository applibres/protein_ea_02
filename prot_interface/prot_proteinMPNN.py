import torch
import torch.nn.functional as F
import numpy as np
from ProteinMPNN.protein_mpnn_utils import ProteinMPNN, parse_PDB

class ProteinMPNNProbMatrix:
    def __init__(self, pdb_path, checkpoint_path="ProteinMPNN/vanilla_model_weights/v_48_030.pt"):
        self.device = torch.device(
            "cuda" if torch.cuda.is_available()
            else "cpu"
        )

        self.AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

        pdb_dict_list = parse_PDB(pdb_path)
        pdb_dict = pdb_dict_list[0]   # tu función siempre devuelve lista
        
        # Detectar cadenas reales presentes
        chains = [
            k.split("_")[-1]
            for k in pdb_dict.keys()
            if k.startswith("seq_chain_")
        ]
        
        # Usar todas las cadenas (IMPORTANTE para interfaz)
        X_all = []
        seq_concat = ""
        
        for chain_id in chains:
            seq = pdb_dict[f"seq_chain_{chain_id}"]
            coords = pdb_dict[f"coords_chain_{chain_id}"]
        
            N  = np.array(coords[f"N_chain_{chain_id}"])
            CA = np.array(coords[f"CA_chain_{chain_id}"])
            C  = np.array(coords[f"C_chain_{chain_id}"])
            O  = np.array(coords[f"O_chain_{chain_id}"])
        
            X_chain = np.stack([N, CA, C, O], axis=1)
        
            X_all.append(X_chain)
            seq_concat += seq
        
        # Concatenar todas las cadenas
        X = np.concatenate(X_all, axis=0)
        
        self.X = torch.tensor(X, dtype=torch.float32).unsqueeze(0).to(self.device)
        
        # mask = 1 para todos los residuos existentes
        self.chain_mask = torch.ones(
            (1, X.shape[0]), dtype=torch.float32
        ).to(self.device)
        
        self.native_seq = seq_concat

        # Load model
        self.model = ProteinMPNN(
            num_letters=21,
            node_features=128,
            edge_features=128,
            hidden_dim=128,
            num_encoder_layers=3,
            num_decoder_layers=3,
            k_neighbors=32,
            dropout=0.1
        ).to(self.device)

        checkpoint = torch.load(checkpoint_path, map_location=self.device)
        self.model.load_state_dict(checkpoint["model_state_dict"])
        self.model.eval()

        self.cache = {}

    def seq_to_tensor(self, seq):
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        return torch.tensor(
            [alphabet.index(a) for a in seq],
            device=self.device
        ).unsqueeze(0)

    @torch.no_grad()
    def _forward(self, seq):
    
        if seq in self.cache:
            return self.cache[seq]
    
        S = self.seq_to_tensor(seq)
    
        B, L = S.shape
        device = self.device
    
        mask = torch.ones((B, L), device=device)
        chain_M = torch.ones((B, L), device=device)
    
        residue_idx = torch.arange(L, device=device).unsqueeze(0)
    
        chain_encoding_all = torch.zeros((B, L), device=device)
    
        randn = torch.randn((B, L), device=device)
    
        log_probs = self.model(
            X=self.X,
            S=S,
            mask=mask,
            chain_M=chain_M,
            residue_idx=residue_idx,
            chain_encoding_all=chain_encoding_all,
            randn=randn
        )
    
        log_probs = log_probs[0, :, :20]
        probs = torch.exp(log_probs)
    
        seq_indices = S[0]
        ll_positions = log_probs[torch.arange(L), seq_indices]
        ll_total = ll_positions.sum().item()
    
        data = {
            "ll": ll_total,
            "ll_positions": ll_positions,
            "prob_matrix": probs
        }
    
        self.cache[seq] = data
        return data

    def get_mpnn_ll(self, seq):
        return self._forward(seq)["ll"]

    def get_probability_matrix(self, seq):
        return self._forward(seq)["prob_matrix"]

    def most_probable_replacement(self, seq, position):
        data = self._forward(seq)
        probs = data["prob_matrix"][position].clone()

        wt = seq[position]
        wt_idx = self.AMINO_ACIDS.index(wt)
        probs[wt_idx] = 0.0

        best_idx = probs.argmax().item()
        return self.AMINO_ACIDS[best_idx], probs[best_idx].item()
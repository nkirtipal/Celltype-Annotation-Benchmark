import os
os.environ["OMP_NUM_THREADS"] = "14"
os.environ["MKL_NUM_THREADS"] = "14"

import argparse, numpy as np, pandas as pd, scanpy as sc, torch
import torch.nn as nn
from pathlib import Path
from sklearn.model_selection import train_test_split

torch.set_num_threads(14)

class ACTINN(nn.Module):
    def __init__(self, n_in, n_out):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_in, 100), nn.ReLU(),
            nn.Linear(100, 50),   nn.ReLU(),
            nn.Linear(50, 25),    nn.ReLU(),
            nn.Linear(25, n_out))
    def forward(self, x):
        return self.net(x)

def preprocess(X):
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    X = X / X.sum(axis=1, keepdims=True) * 1e4
    return np.log2(X + 1).astype(np.float32)

def select_genes(Xl):
    expr = Xl.sum(axis=0)
    lo, hi = np.percentile(expr, [1, 99])
    m = (expr >= lo) & (expr <= hi)
    Xl, idx = Xl[:, m], np.where(m)[0]
    mu = Xl.mean(axis=0); keep = mu > 0
    Xl, idx = Xl[:, keep], idx[keep]
    cv = Xl.std(axis=0) / Xl.mean(axis=0)
    lo, hi = np.percentile(cv, [1, 99])
    k2 = (cv >= lo) & (cv <= hi)
    return idx[k2]

def main(seed):
    ROOT = Path("/media/nikhil/My_Book/celltype_bench")
    ad = sc.read_h5ad(ROOT / "data/processed/pbmc3k.h5ad")

    y_str = ad.obs["cell_type"].astype(str).values
    classes = sorted(set(y_str))
    y = np.array([classes.index(c) for c in y_str])

    tr, te = train_test_split(np.arange(ad.n_obs), test_size=0.3,
                              stratify=y, random_state=seed)

    Xl = preprocess(ad.X)
    gidx = select_genes(Xl[tr])          # gene selection on TRAIN only
    Xtr, Xte = Xl[tr][:, gidx], Xl[te][:, gidx]

    torch.manual_seed(seed); np.random.seed(seed)
    model = ACTINN(len(gidx), len(classes))
    opt = torch.optim.Adam(model.parameters(), lr=1e-4, weight_decay=5e-3)
    sched = torch.optim.lr_scheduler.ExponentialLR(opt, gamma=0.95)
    lossf = nn.CrossEntropyLoss()

    Xt = torch.tensor(Xtr); yt = torch.tensor(y[tr])
    for ep in range(50):
        perm = torch.randperm(len(Xt))
        for i in range(0, len(Xt), 128):
            b = perm[i:i+128]
            opt.zero_grad()
            l = lossf(model(Xt[b]), yt[b]); l.backward(); opt.step()
        sched.step()
        if (ep+1) % 10 == 0:
            print(f"  epoch {ep+1} loss {l.item():.4f}")

    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(Xte))
        prob = torch.softmax(logits, 1).numpy()

    pd.DataFrame({
        "barcode": ad.obs_names[te],
        "label": [classes[i] for i in prob.argmax(1)],
        "confidence": prob.max(1),
        "true_label": y_str[te],
    }).to_csv(ROOT / f"results/predictions_actinn_seed{seed}_pbmc3k.csv", index=False)

    acc = (prob.argmax(1) == y[te]).mean()
    print(f"seed {seed}: n_genes={len(gidx)} test_acc={acc:.4f}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--seed", type=int, default=0)
    main(p.parse_args().seed)

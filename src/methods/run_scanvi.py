import os
os.environ["OMP_NUM_THREADS"] = "14"
os.environ["MKL_NUM_THREADS"] = "14"

import argparse, numpy as np, pandas as pd, scanpy as sc, scvi, torch
from pathlib import Path
from sklearn.model_selection import train_test_split

torch.set_num_threads(14)

def main(seed):
    ROOT = Path("/media/nikhil/My_Book/celltype_bench")
    ad = sc.read_h5ad(ROOT / "data/processed/pbmc3k.h5ad")
    scvi.settings.seed = seed

    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, subset=True,
                                layer="counts", flavor="seurat_v3")
    ad.X = ad.layers["counts"].copy()

    y = ad.obs["cell_type"].astype(str).values
    tr, te = train_test_split(np.arange(ad.n_obs), test_size=0.3,
                              stratify=y, random_state=seed)

    ad.obs["labels_semi"] = "Unknown"
    ad.obs.iloc[tr, ad.obs.columns.get_loc("labels_semi")] = y[tr]

    scvi.model.SCVI.setup_anndata(ad, layer="counts", labels_key="labels_semi")
    vae = scvi.model.SCVI(ad, n_latent=30)
    vae.train(max_epochs=200, accelerator="cpu", early_stopping=True)

    lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="Unknown")
    lvae.train(max_epochs=100, accelerator="cpu", early_stopping=True)

    pred = lvae.predict(ad)
    prob = lvae.predict(ad, soft=True).max(axis=1).values

    out = pd.DataFrame({
        "barcode": ad.obs_names, "label": pred,
        "confidence": prob, "true_label": y,
        "split": np.where(np.isin(np.arange(ad.n_obs), tr), "train", "test"),
    })
    out.to_csv(ROOT / f"results/predictions_scanvi_seed{seed}_pbmc3k.csv", index=False)

    t = out[out.split == "test"]
    print(f"seed {seed}: test_acc={(t.label == t.true_label).mean():.4f}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--seed", type=int, default=0)
    main(p.parse_args().seed)

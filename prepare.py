import os
import cellxgene_census


DATA_DIR = "data"
OUT_PATH = os.path.join(DATA_DIR, "lims_lung_celltype_demo.h5ad")


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    if os.path.exists(OUT_PATH):
        print(f"Dataset already exists: {OUT_PATH}")
        return

    keep_cell_types = [
        "T cell",
        "B cell",
        "natural killer cell",
        "macrophage",
        "monocyte",
        "dendritic cell",
        "epithelial cell",
        "endothelial cell",
        "fibroblast",
    ]

    with cellxgene_census.open_soma(census_version="stable") as census:
        obs = cellxgene_census.get_obs(
            census,
            organism="Homo sapiens",
            value_filter=(
                "tissue_general == 'lung' "
                "and is_primary_data == True "
                "and disease in ['normal', 'COVID-19']"
            ),
            column_names=[
                "soma_joinid",
                "cell_type",
                "disease",
                "tissue",
                "tissue_general",
                "assay",
                "dataset_id",
                "donor_id",
            ],
        )

        obs = obs[obs["cell_type"].isin(keep_cell_types)].copy()

        # Keep the demo small and balanced enough for a classroom run.
        obs_small = (
            obs.groupby("cell_type", group_keys=False)
            .apply(lambda x: x.sample(min(len(x), 500), random_state=1))
            .reset_index(drop=True)
        )

        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            obs_coords=obs_small["soma_joinid"].to_numpy(),
            column_names={
                "obs": [
                    "cell_type",
                    "disease",
                    "tissue",
                    "tissue_general",
                    "assay",
                    "dataset_id",
                    "donor_id",
                ]
            },
        )

    adata.write_h5ad(OUT_PATH)
    print(adata)
    print(f"Wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
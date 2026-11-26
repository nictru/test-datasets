# scdownstream atlas extension fixture

The merged H5AD and matching scVI and Symphony models come from the default
pipeline nf-test run at scdownstream commit `57c37bf63186c5686757b1b136752dd297aecae2`.
Run `nf-test test tests/default.nf.test` with a container profile on that revision.
Copy `09_finalized/merged.h5ad`, `03_combine/integrate/scvi/scvi_model/model.pt`,
and `03_combine/integrate/symphony/symphony_reference.h5ad` from the test output.
Remove the `uns/gene_presence_complete` entry from the H5AD with h5py. The
per-sample presence table comes from the de novo pipeline run.

FoldBounds — exact GNM and thioredoxin blind experiment

STAGE 1 (NEW, no target coordinates):
Download only the file fb_thioredoxin_stage1.py to your Windows Downloads folder.
In PyMOL 3.1.8:
    run C:/Users/amarena/Downloads/fb_thioredoxin_stage1.py
    fb_thio_freeze

PyMOL loads ONLY X-ray 2TRX chain A and freezes all model arrays.
Send back the entire generated FoldBounds_THIO_STAGE1_20261105_SEND_BACK.zip.
DO NOT run Stage 2 or fetch the NMR target 1XOA until Stage1 SHA is audited.

STAGE 2 (CODE PREPARED, DO NOT RUN YET):
    run C:/Users/amarena/Downloads/fb_thioredoxin_stage2.py
    fb_thio_evaluate <EXTERNALLY VERIFIED STAGE1 ZIP SHA256>

Stage2 MUST reject missing/incorrect 64-char Stage1 SHA prior to NMR access.
Predeclared primary model: generalized ANM, 50/50 isotropic/directional.
Criterion: rho>=0.50 AND +0.15 over strongest of SIX old baseline models.
The physical model is not novel: see Zheng 2008 G-ANM literature.

Independent reproducibility archives:
  gnm_rnase_reference_certificate.json (2,614 reference-only GNM interval certificates)
  gnm_structural_robustness.json (all optional contact graph states bounded)
  scripts/gnm_independent_replay.py and gnm_independent_box_replay.py
  research/ExactGNM_Topology_UQ_20261105.md (full audit)

Never treat the old RNase A failure as a success or mistake. No 1XOA targets
were used in choosing model settings. NMR model spread is NOT physical dynamics.
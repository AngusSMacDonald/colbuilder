Colbuilder’s coarse-grained build path is:

1. Geometry/mixing writes per-helix caps files (_*.caps.pdb_) and connectivity (_connect_from_colbuilder.txt_) via Connect.write\_connect (_connect.py_ (line 137), _connect.py_ (line 156)).
2. The mix helper reconstructs a System by assigning every ID on a line the same connect list (_topology_from_mix.py_ (line 44), _topology_from_mix.py_ (line 51), _topology_from_mix.py_ (line 62)).
3. Topology then runs build\_martini3 (_main_topology.py_ (line 202)), which loops over each model and each connect\_id (_martini.py_ (line 722), _martini.py_ (line 744)), reads _<connect_id>.caps.pdb_ (_martini.py_ (line 187)), writes _<model_id>.<connect_id>.CG.pdb_ (_martini.py_ (line 762)), merges to _*.merge.pdb_ (_martini.py_ (line 821)), then emits _col\_\<cnt\_model\>.itp_ (_itp.py_ (line 643)).

Why your overlaps happen:

- If two models have the same connection set (example [8, 38]), both models run the same connect\_id inputs.
- The coordinates come from connect\_id caps files, not from current model\_id (martini.py (line 187)), so both builds produce the same CG coordinates.
- They are written under different molecule indices (col\_8, col\_38, etc.), so you see bead-for-bead overlapping duplicates with different names.
This is not just the helper script; the original topology path also calls the same build\_martini3 loop, so the duplication mechanism is shared.

I also verified this behavior in repo outputs: in runs/test/.tmp/topology\_gen, multiple _*.merge.pdb_ files are identical pairs (e.g., 16.merge.pdb = 17.merge.pdb, 6.merge.pdb = 7.merge.pdb).



### This is a core code issue, needing changes in martini.py ###

Implemented the core fix in martini.py to deduplicate duplicate coarse-grained builds during Martini topology generation.

#What I changed:#

1. Added connection normalization helper in build\_martini3 at _martini.py_ (line 723).
 - New \_normalize\_connect\_ids(...) converts connection IDs to deterministic ints, removes duplicates, sorts them, and falls back to self-connection if needed.
 - This avoids group-key mismatches from float/string ID formatting and repeated IDs.
2. Added group-level deduplication state at martini.py (line 721).
 - New processed\_groups: Dict[Tuple[str, Tuple[int, ...]], float].
3. Replaced old direct connect handling with normalized IDs at martini.py (line 753).
 - Old logic used raw model.connect or [model\_id].
 - New logic uses \_normalize\_connect\_ids(...) consistently.
4. Added skip logic for duplicate groups at martini.py (line 762).
 - Group key is (model.type, tuple(connect_ids)).
 - If already processed, model is skipped and marked in model_status as duplicate_group_of_<representative>.

Previously, multiple models sharing the same connection component were all re-martinized, producing identical merged structures under different col_N labels.
Now, each unique (type + connection set) is built once, preventing overlapping duplicate molecules like the col_8/col_38 pattern you described.

### THE REAL ISSUE ###
There are only Martini3 parameterisations done for PYD and HLKNL crosslinks (which would be so great if they mentioned that somewhere...). Therefore I need to parameterise DPD, LKNL, deHHLNL and deHLNL myself.

DPD is hydrolysed version of PYD - ie OH group on backbone is removed
LKNL is hydrolysed version of HLKNL

deHHLNL is a de-oxidised (?) version of HLKNL, and deHLNL is that for LKNL (_double-bonded oxygen removed, nitrogen on backbone double-bonded to instead_)

I have made temp files for DPD and LKNL, placeholders for now but topology will build
I need to do this for deHHLNL and deHLNL still.

ALL NEW PARAMETERS STILL NEED TO BE ADDED. CURRENTLY PLACEHOLDERS (_WRONG_) OR MISSING!

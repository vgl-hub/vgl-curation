# CHANGELOG


## v1.1 - 2026-01-21

### Bug Fixes

- **split_agp.py**: Fixed component numbering for unloc sequences. Component numbers (column 4) now correctly reset to 1 for each new unloc object instead of inheriting numbering from parent scaffold.
- **split_agp.py**: Fixed object coordinate calculation for unloc sequences. Object end coordinates (column 3) now correctly calculate length from component coordinates `(comp_end - comp_start + 1)` instead of incorrectly using the raw component end coordinate. This resolves AGP validation errors: "object coordinates and component coordinates do not have the same length".

### Technical Details

Fixed in `unloc()` function (lines 269-273):
- Added: `agp_df.loc[index,'#_scaffs']=1` to reset component numbering
- Changed: `agp_df.loc[index,'chr_end']` calculation from `scaff_end` to `int(scaff_end) - int(scaff_start) + 1`


## v1.0 - 2025-08-06

- Initial Release




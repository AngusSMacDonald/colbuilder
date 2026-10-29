#!/usr/bin/env bash
# Build one ColBuilder system and check it with GROMACS: vacuum EM -> solvate ->
# genion -> solvated EM. Minimisations are shortened and the box is fitted to the
# fibril (-d 1.5). The python checks in this directory then run on the topology.
#
# usage: run_topology_check.sh TEST_DIR
#   TEST_DIR  directory holding config.yaml; all outputs are written there
#
# Optional environment:
#   GMXRC          GROMACS environment script to source (otherwise gmx must be on PATH)
#   GMX_THREADS    OpenMP threads for mdrun (default 8)
#   VALIDATION_LOCK  file to flock around the GROMACS steps, so concurrent checks do
#                  not overlap; if unset, no lock is taken
#
# MODELLER and colbuilder must already be set up as in the README, so that
# `colbuilder` is on PATH in the activated environment.
#
# Writes TEST_DIR/STATUS with one "STEP OK|FAIL detail" line per stage.
set -uo pipefail

[ $# -eq 1 ] || { echo "usage: run_topology_check.sh TEST_DIR" >&2; exit 2; }
T=$(cd "$1" && pwd) || { echo "no such directory: $1" >&2; exit 2; }
H=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GMX_THREADS=${GMX_THREADS:-8}
if [ -n "${GMXRC:-}" ]; then set +u; source "$GMXRC"; set -u; fi

ST="$T/STATUS"; cd "$T"; : > "$ST"
ok()   { echo "$1 OK${2:+ $2}" >> "$ST"; }
bad()  { echo "$1 FAIL ${2:-}" >> "$ST"; checks_failed=1; }
fail() { echo "$1 FAIL ${2:-}" >> "$ST"; echo "DONE FAIL" >> "$ST"; exit 1; }
em_summary() {  # final energy / max force from an EM log
    awk '/^Potential Energy/{e=$4} /^Maximum force/{f=$4} END{printf "Epot=%s Fmax=%s", e, f}' "$1"
}
checks_failed=0

start=$SECONDS
( cd "$T" && colbuilder --config_file config.yaml ) > build.log 2>&1 \
    || fail BUILD "exit $? after $((SECONDS-start))s; see build.log"
if grep -q "Unexpected error" build.log; then fail BUILD "'Unexpected error' in build.log"; fi
top=$(find "$T" -maxdepth 2 -name 'collagen_fibril_*.top' -not -path '*/.tmp/*' | head -1)
[ -n "$top" ] || fail BUILD "no .top produced; see build.log"
ok BUILD "$((SECONDS-start))s top=${top#$T/}"

# Everything below runs in the .top directory so relative includes resolve.
D=$(dirname "$top"); cd "$D"; stem=$(basename "$top" .top)
[ -f "$stem.gro" ] || fail BUILD "no $stem.gro next to $stem.top"

if [ -n "${VALIDATION_LOCK:-}" ]; then exec 9>"$VALIDATION_LOCK"; flock 9; fi
sed 's/^\s*nsteps.*/nsteps = 500/' "$H/mdps/minim_vacuum.mdp" > em_vac.mdp
sed 's/^\s*nsteps.*/nsteps = 200/' "$H/mdps/minim_solv.mdp" > em_solv.mdp

gmx editconf -f "$stem.gro" -o boxed.gro -c -d 1.5 > editconf.log 2>&1 || fail EDITCONF "see editconf.log"
ok EDITCONF "$(tail -1 boxed.gro | xargs)"

gmx grompp -f em_vac.mdp -c boxed.gro -p "$stem.top" -o em_vac.tpr -maxwarn 1 > grompp_vac.log 2>&1 \
    || fail GROMPP_VAC "see $D/grompp_vac.log"
ok GROMPP_VAC "$(grep -c '^WARNING' grompp_vac.log) warnings, $(grep -c '^NOTE' grompp_vac.log) notes"
gmx mdrun -deffnm em_vac -ntmpi 1 -ntomp "$GMX_THREADS" > mdrun_vac.out 2>&1 || fail EM_VAC "see $D/mdrun_vac.out"
ok EM_VAC "$(em_summary em_vac.log)"

cp "$stem.top" solv.top
if grep -q 'tip3p.itp' "$stem.top"; then water=spc216.gro; else water=tip4p.gro; fi
gmx solvate -cp em_vac.gro -cs "$water" -p solv.top -o solv.gro > solvate.log 2>&1 || fail SOLVATE "see solvate.log"
ok SOLVATE "$(sed -n 2p solv.gro | xargs) atoms ($water)"

gmx grompp -f "$H/mdps/ions.mdp" -c solv.gro -p solv.top -o ions.tpr > grompp_ions.log 2>&1 \
    || fail GROMPP_IONS "see $D/grompp_ions.log"
cp solv.top ions.top
printf 'SOL\n' | gmx genion -s ions.tpr -o ions.gro -p ions.top -conc 0.15 -neutral > genion.log 2>&1 \
    || fail GENION "see genion.log"
ok GENION "$(grep -E 'Will try to add' genion.log | xargs)"

gmx grompp -f em_solv.mdp -c ions.gro -p ions.top -o em_solv.tpr -maxwarn 1 > grompp_solv.log 2>&1 \
    || fail GROMPP_SOLV "see $D/grompp_solv.log"
ok GROMPP_SOLV "$(grep -c '^WARNING' grompp_solv.log) warnings, $(grep -c '^NOTE' grompp_solv.log) notes"
gmx mdrun -deffnm em_solv -ntmpi 1 -ntomp "$GMX_THREADS" > mdrun_solv.out 2>&1 || fail EM_SOLV "see $D/mdrun_solv.out"
ok EM_SOLV "$(em_summary em_solv.log)"

rm -f solv.gro ions.tpr em_vac.trr em_solv.trr   # large intermediates
exec 9>&-   # GROMACS is finished; release the lock before the python checks

# Python checks. A FAIL is recorded but the remaining checks still run.
unbonded_summary() { grep -o 'close_pairs=.*' unbonded.txt | head -1 || tail -1 unbonded.txt; }
if python3 "$H/unbonded_xlinks.py" "$stem.top" "$stem.gro" . > unbonded.txt 2>&1; then
    ok UNBONDED "$(unbonded_summary)"
else
    bad UNBONDED "$(unbonded_summary)"
fi
if python3 "$H/xlink_terms.py" . > xlink_terms.txt 2>&1; then
    ok XLINK_TERMS
else
    bad XLINK_TERMS "see $D/xlink_terms.txt"
fi
if python3 "$H/xlink_geometry.py" "$stem.top" . "$stem.gro" em_vac.gro > xlink_geometry.txt 2>&1; then
    ok XLINK_GEOMETRY
else
    bad XLINK_GEOMETRY "see $D/xlink_geometry.txt"
fi

if [ "$checks_failed" -eq 0 ]; then
    echo "DONE OK $((SECONDS-start))s" >> "$ST"
else
    echo "DONE FAIL $((SECONDS-start))s" >> "$ST"
    exit 1
fi

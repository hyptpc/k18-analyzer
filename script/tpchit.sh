#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."
git checkout tpchit/v1

# =================================================
EXAMPLE_USER_CC="example/UserTPCHit.cc"  
USR_DIR="usr"
USR_USER_CC="${USR_DIR}/UserTPCHit.cc"

ASSIGNMENT_FILE="runmanager/runlist/assignment_summary.txt"
YML_OUTDIR="runmanager/runlist"


MAKEFILE_ORG="Makefile.org"
MAKEFILE="Makefile"

QUEUE="l"
QMERGE="l"
UNIT=1000
NPROC=18

BUFF="/ghi/fs02/orig_root_fs02/had/sks/Users/E72_TPCHit"
BIN="./bin/TPCHit"
CONF="param/conf/analyzer_e72_tpc_0.conf"
DATA="/hsm/had/sks/E72/JPARC2025Nov/e72_2025nov"
ROOT_BASE="/hsm/had/sks/E72/JPARC2025Nov/rootfile/tpchit/v1"
# =================================================

usage () {
  echo "Usage: $0 <your_name> [--make-j N] [--session NAME] [--no-run]"
  exit 1
}

[[ $# -lt 1 ]] && usage

NAME="$1"; shift
MAKEJ=20
SESSION="tpchit_${NAME}"
DO_RUN=1


while [[ $# -gt 0 ]]; do
  case "$1" in
    --make-j) MAKEJ="$2"; shift 2;;
    --session) SESSION="$2"; shift 2;;
    --no-run) DO_RUN=0; shift;;
    *) usage;;
  esac
done

if [[ -z "${TMUX:-}" ]]; then
  echo "[INFO] Creating tmux session: ${SESSION}"
  tmux new-session -d -s "${SESSION}"
  tmux send-keys -t "${SESSION}" \
    "$0 ${NAME} --make-j ${MAKEJ}" C-m
  tmux attach -t "${SESSION}"
  exit 0
fi

YML_PATH="runmanager/runlist/${NAME}.yml"
echo "[INFO] tmux session: ${SESSION}"
echo "[INFO] WORKDIR = $(pwd)"

# ---------- UserTPCHit.cc ----------
mkdir -p "${USR_DIR}"
cp -f "${EXAMPLE_USER_CC}" "${USR_USER_CC}"
echo "[OK] UserTPCHit.cc prepared"

# ---------- Makefile ----------
if [[ -f "${MAKEFILE_ORG}" ]]; then
  cp -f "${MAKEFILE_ORG}" "${MAKEFILE}"
  echo "[OK] Makefile.org -> Makefile"
fi

# ---------- Build ----------
echo "[INFO] make -j ${MAKEJ}"
make -j "${MAKEJ}"
echo "[OK] build done"

# ---------- YML ----------
mkdir -p "${YML_OUTDIR}"

python3 - "$NAME" "$(pwd)" "$YML_PATH" <<'PY'
import re, sys
from pathlib import Path

name = sys.argv[1]
workdir = sys.argv[2]
yml_path = Path(sys.argv[3])

assign = "runmanager/runlist/assignment_summary.txt"
outdir = Path("runmanager/runlist")

QUEUE="l"; QMERGE="l"
UNIT=1000; NPROC=18
BUFF="/ghi/fs02/orig_root_fs02/had/sks/Users/E72_TPCHit"
BIN="./bin/TPCHit"
CONF="param/conf/analyzer_e72_tpc_0.conf"
DATA="/hsm/had/sks/E72/JPARC2025Nov/e72_2025nov"
ROOT_BASE="/hsm/had/sks/E72/JPARC2025Nov/rootfile/tpchit/v1"

runs = None
with open(assign) as f:
    for line in f:
        m = re.match(rf"{name}\s*\(\d+\):\s*(.*)", line.strip())
        if m:
            tail = m.group(1)
            runs = [int(x) for x in tail.split(",")] if tail else []
            break

if runs is None:
    sys.exit("Name not found in assignment_summary.txt")

root = Path(ROOT_BASE) / name
fig  = root / "fig"
root.mkdir(parents=True, exist_ok=True)
fig.mkdir(exist_ok=True)

run_lines = "".join(f"  {r}:\n" for r in runs)

yml = f"""WORKDIR: {workdir}

DEFAULT:
  queue:  {QUEUE}
  qmerge: {QMERGE}
  unit:   {UNIT}
  nproc:  {NPROC}
  buff:   {BUFF}
  bin:    {BIN}
  conf:   {CONF}
  data:   {DATA}
  root:   {root}
  fig:    {fig}

RUN:
{run_lines}
"""

yml_path.parent.mkdir(parents=True, exist_ok=True)
yml_path.write_text(yml)
print(yml_path)
PY

echo "[OK] YML created: ${YML_PATH}"

# ---------- Prepare output directory ----------
ROOT_DIR="${ROOT_BASE}/${NAME}"
FIG_DIR="${ROOT_DIR}/fig"

mkdir -p "${ROOT_DIR}"
mkdir -p "${FIG_DIR}"

echo "[OK] Output directory prepared:"
echo "     ${ROOT_DIR}"

# ---------- run.py ----------
if [[ "${DO_RUN}" -eq 0 ]]; then
  echo "[INFO] --no-run : stop here"
  exit 0
fi


echo "[INFO] Running:"
echo "  ./runmanager/run.py ${YML_PATH}"

set -o pipefail
./runmanager/run.py "${YML_PATH}" 2>&1

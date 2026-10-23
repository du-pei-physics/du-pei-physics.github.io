#!/bin/zsh
# Benchmark v3 (19 Oct 2026), after Codex's publication review: every hit is verified from the generator and
# hashed (distinct rowspaces across pools and runs), search and verification CPU are recorded separately, pools are
# reproducible (RNG state per pool in epochs.jsonl), and every run's seed, argv and source hashes are in runs.jsonl.
# Stopping: touch persistent/bench3_stop.flag stops every lane within seconds (workers poll it); a lane's own
# persistent/bench3_kX/stop.flag stops that configuration and is honoured, never deleted, by this wrapper.
# A failing child is logged in persistent/bench3_failures.log; after three consecutive failures the lane exits.
# Lane G (19 Oct): 20 minutes of Q8-equivariant chains at k = 8 (bench_q8.py, modular pre-screen) before each k = 21 slot;
# the k = 8 slot is paused since 19 Oct 07:30 by bench3_q8_k8/stop.flag (4,884 codes; workers moved to few-hit k).
# Lanes A, B since 19 Oct 07:30: k = 9 (chains of cells 12, 6, period 3, three boundary spins) and k = 16 (Q8 chains);
# since 19 Oct 12:00: A = third k = 18 lane, B = k = 16 with the batched block engine (estimates still from few codes).
# usage: bench3_run.sh LANE   (LANE = A..G)
cd "$(dirname "$0")/mixed_rings" || exit 1
STOP="$(cd .. && pwd)/bench3_stop.flag"
FAILLOG="$(cd .. && pwd)/bench3_failures.log"
P=../../../.venv/bin/python
COMMON=(--exclude-bases A0 --s-refl 2 --pool 40 --refresh 200000 --hits-max 100 --global-stop "$STOP")
fails=0
run() { # name k widths periods r fast metric2 minutes
  [ -f "$STOP" ] && exit 0
  if [ -f "../$1/stop.flag" ]; then sleep 60; return 0; fi   # honour a lane stop without spinning
  mkdir -p "../$1"
  local seed=$(( (RANDOM << 15) | RANDOM ))
  nice -n 19 $P bench_worker.py --name $1 --k $2 --widths $3 --periods $4 --r $5 --fast $6 --metric2 $7 --max-minutes $8 --seed $seed $COMMON >> "../$1/log.txt" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') lane $LANE run $1 exit $rc (seed $seed)" >> "$FAILLOG"
    fails=$((fails + 1)); [ $fails -ge 3 ] && exit 1; sleep 60
  else
    fails=0
  fi
}
run_q8() { # name k minutes  (Q8-equivariant chains, bench_q8.py; added 19 Oct for the k=8 yield)
  [ -f "$STOP" ] && exit 0
  if [ -f "../$1/stop.flag" ]; then sleep 60; return 0; fi
  mkdir -p "../$1"
  local seed=$(( (RANDOM << 15) | RANDOM ))
  nice -n 19 $P bench_q8.py --name $1 --k $2 --max-minutes $3 --seed $seed --hits-max 100 --global-stop "$STOP" >> "../$1/log.txt" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') lane $LANE run $1 exit $rc (seed $seed)" >> "$FAILLOG"
    fails=$((fails + 1)); [ $fails -ge 3 ] && exit 1; sleep 60
  else
    fails=0
  fi
}
run_q8b() { # name k minutes  (Q8 chains, batched block engine bench_q8_blocks.py; added 19 Oct 12:00)
  [ -f "$STOP" ] && exit 0
  if [ -f "../$1/stop.flag" ]; then sleep 60; return 0; fi
  mkdir -p "../$1"
  local seed=$(( (RANDOM << 15) | RANDOM ))
  nice -n 19 $P bench_q8_blocks.py --name $1 --k $2 --max-minutes $3 --seed $seed --hits-max 100 --global-stop "$STOP" >> "../$1/log.txt" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') lane $LANE run $1 exit $rc (seed $seed)" >> "$FAILLOG"
    fails=$((fails + 1)); [ $fails -ge 3 ] && exit 1; sleep 60
  else
    fails=0
  fi
}
LANE=$1
while true; do
  [ -f "$STOP" ] && exit 0
  case $LANE in
    A) run bench3_k18c 18 12,12,12 3 4 1 1 360 ;;   # 19 Oct 12:00: k = 9 done (391 codes, estimate settled); third k = 18 lane (estimate from one code)
    B) run_q8b bench3_q8b_k16 16 360 ;;      # 19 Oct 12:00: k = 16 with the batched block engine (bench_q8_blocks.py, ~2-3x faster)
    C) run bench3_k21c 21 12,12,12,6 3 4 1 1 360 ;;   # 19 Oct 15:30: k = 15 done (46 codes, estimate within ~15%); third k = 21 lane
    D) run bench3_k18 18 12,12,12 3 4 1 1 360 ;;
    E) run bench3_k21 21 12,12,12,6 3 4 1 1 360 ;;
    F) run bench3_k18m1 18 12,12,12 3 4 1 0 360 ;;
    G) run_q8 bench3_q8_k8 8 20; run bench3_k21b 21 12,12,12,6 3 4 1 1 360 ;;
    *) echo "unknown lane $LANE"; exit 1 ;;
  esac
done

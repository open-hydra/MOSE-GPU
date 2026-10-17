#!/bin/bash
# mose_bind.sh -- per-rank CPU placement for MOSE_GPU (speed-up 30/09). Usage: mpirun -n N --bind-to none mose_bind.sh ./bin/MOSE [args]
# MOSE gives node-local rank k the CUDA device k mod ndev (Mod_MPI.f90, mpi_bind_device). CUDA device k is the k-th entry of
# CUDA_VISIBLE_DEVICES (set by Slurm; an index or a UUID), read in PCI order: the script exports CUDA_DEVICE_ORDER=PCI_BUS_ID so that
# CUDA's order and nvidia-smi's are the same by construction. The rank then runs on the PHYSICAL cores that are both allowed to the job
# (its cgroup cpuset) and local to its GPU (/sys/bus/pci/devices/<bus>/local_cpulist); the ranks whose GPUs share a local set split its
# cores into disjoint equal blocks (all hardware threads of a core go together). Fewer than MOSE_BIND_MIN cores (default 4) for a rank, or anything unreadable: a warning, and the rank runs
# unbound, exactly as with --bind-to none. MOSE_BIND=0 disables it. Results do not depend on the placement.
# Launchers: OpenMPI mpirun (OMPI_COMM_WORLD_LOCAL_RANK/SIZE), Intel MPI and MPICH hydra (MPI_LOCALRANKID/MPI_LOCALNRANKS), srun with
# one task per rank (SLURM_LOCALID + SLURM_STEP_TASKS_PER_NODE); with anything else it runs unbound and says so.
export CUDA_DEVICE_ORDER=PCI_BUS_ID
[ "${MOSE_BIND:-1}" = 0 ] && exec "$@"
rk=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${SLURM_PROCID:-?}}}
unbound() { echo "[mose_bind] rank $rk: $1 -> unbound (as --bind-to none)" >&2; exec "${ARGS[@]}"; }
ARGS=("$@")
expand() { local out=() a b p; IFS=',' read -ra parts <<< "$1"; for p in "${parts[@]}"; do
             if [[ $p == *-* ]]; then a=${p%-*}; b=${p#*-}; for ((c=a; c<=b; c++)); do out+=($c); done; elif [ -n "$p" ]; then out+=($p); fi
           done; echo "${out[@]}"; }
# the rank's place on this node, by launcher: OpenMPI mpirun; Intel MPI / MPICH hydra; srun (one task per rank). Anything else:
# unbound, because a wrong local rank would put every rank on the cores of GPU 0 (silent, results unchanged, much slower).
if [ -n "${OMPI_COMM_WORLD_LOCAL_RANK:-}" ]; then lr=$OMPI_COMM_WORLD_LOCAL_RANK; ln=${OMPI_COMM_WORLD_LOCAL_SIZE:-}
elif [ -n "${MPI_LOCALRANKID:-}" ]; then lr=$MPI_LOCALRANKID; ln=${MPI_LOCALNRANKS:-}
elif [ -n "${SLURM_LOCALID:-}" ] && [ -n "${SLURM_STEP_TASKS_PER_NODE:-}" ]; then
  lr=$SLURM_LOCALID; ln=${SLURM_STEP_TASKS_PER_NODE%%,*}; ln=${ln%%(*}
else unbound "unknown launcher: no local rank in the environment"; fi
[[ $lr =~ ^[0-9]+$ && $ln =~ ^[0-9]+$ ]] || unbound "local rank '$lr' / local size '$ln' not numbers"
[ "$lr" -lt "$ln" ] || unbound "local rank $lr >= local size $ln"
allowed=" $(expand "$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)") "
# the visible devices, in CUDA order (PCI order: CUDA_DEVICE_ORDER above)
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then IFS=',' read -ra vis <<< "$CUDA_VISIBLE_DEVICES"
else mapfile -t vis < <(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null); fi
ngpu=${#vis[@]}; [ "$ngpu" -gt 0 ] || unbound "no GPU visible"
declare -A busof
while IFS=', ' read -r idx bus uuid; do busof[$idx]=$bus; busof[$uuid]=$bus; done < <(nvidia-smi --query-gpu=index,pci.bus_id,uuid --format=csv,noheader 2>/dev/null)
local_of() {   # local cpulist of the k-th visible device
  local id=${vis[$1]} bus sys; bus=${busof[$id]}; [ -n "$bus" ] || return 1
  sys=$(echo "$bus" | tr 'A-F' 'a-f' | sed 's/^0000\(....:\)/\1/')
  cat /sys/bus/pci/devices/$sys/local_cpulist 2>/dev/null; }
myg=$(( lr % ngpu )); mine=$(local_of $myg) || true
[ -n "$mine" ] || unbound "the local CPU list of visible GPU $myg (${vis[$myg]}) is not readable"
share=(); for ((r=0; r<ln; r++)); do [ "$(local_of $(( r % ngpu )))" = "$mine" ] && share+=($r); done
pos=0; for i in "${!share[@]}"; do [ "${share[$i]}" = "$lr" ] && pos=$i; done
# physical cores: CPU -> core id (lscpu); keep a core only if all its CPUs we may use are allowed and GPU-local
declare -A coreof; cores=()
while IFS=',' read -r cpu core sock; do [[ $cpu =~ ^[0-9]+$ ]] && coreof[$cpu]="$sock:$core"; done < <(lscpu -p=CPU,CORE,SOCKET 2>/dev/null)
declare -A cpus_of_core
for c in $(expand "$mine"); do [[ $allowed == *" $c "* ]] || continue; k=${coreof[$c]:-c$c}
  [ -z "${cpus_of_core[$k]}" ] && cores+=("$k"); cpus_of_core[$k]="${cpus_of_core[$k]:+${cpus_of_core[$k]},}$c"; done
per=$(( ${#cores[@]} / ${#share[@]} ))
[ "$per" -ge "${MOSE_BIND_MIN:-4}" ] || unbound "${#cores[@]} allowed physical cores local to GPU $myg for ${#share[@]} rank(s): fewer than ${MOSE_BIND_MIN:-4} each"
mycpus=""; for k in "${cores[@]:$(( pos * per )):$per}"; do mycpus="${mycpus:+$mycpus,}${cpus_of_core[$k]}"; done
echo "[mose_bind] rank $rk (local $lr) -> CUDA device $myg (${vis[$myg]}, bus ${busof[${vis[$myg]}]}), $per cores: cpus $mycpus" >&2
exec taskset -c "$mycpus" "${ARGS[@]}"

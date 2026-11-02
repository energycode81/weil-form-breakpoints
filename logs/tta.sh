# time-to-accuracy: Fourier (weil_step1.py) vs breakpoint (weil_breakpoints.py)
# one line per (lambda, basis, size) in tta.csv: lam,basis,size,dim,seconds,logfile

echo "lam,basis,size,dim,seconds,log" > tta.csv

run () {   # run <lam> <basis> <size> <dim> <cmd...>
  local lam=$1 basis=$2 size=$3 dim=$4; shift 4
  local log="tta_${basis}_lam${lam}_${size}.log"
  local t0=$(date +%s.%N)
  "$@" > "$log" 2>&1
  local t1=$(date +%s.%N)
  echo "$lam,$basis,$size,$dim,$(echo "$t1 - $t0" | bc),$log" >> tta.csv
}

# ---- lambda = 2  (dps 50) ----
for N in 10 15 20 30 40 60 80 100; do
  run 2 fourier $N $((2*N+1)) python3 weil_step1.py --lam 2 --N $N --dps 50 --no-prolate --no-spectrum
done
for d in 4 6 8 10 12 16 20; do
  run 2 breakpoint $d "" python3 weil_breakpoints.py --lam 2 --d $d --dps 50 --ygl 24
done

# ---- lambda = 2.5  (dps 60) ----
for N in 10 15 20 30 40 60 80 100; do
  run 2.5 fourier $N $((2*N+1)) python3 weil_step1.py --lam 2.5 --N $N --dps 60 --no-prolate --no-spectrum
done
for d in 4 6 8 10 12 16 20; do
  run 2.5 breakpoint $d "" python3 weil_breakpoints.py --lam 2.5 --d $d --dps 60 --ygl 24
done

# ---- lambda = 3  (dps 90) -- breakpoint side is slow here ----
for N in 10 20 30 40 60 80 100; do
  run 3 fourier $N $((2*N+1)) python3 weil_step1.py --lam 3 --N $N --dps 90 --no-prolate --no-spectrum
done
for d in 4 6 8 10 12 16; do
  run 3 breakpoint $d "" python3 weil_breakpoints.py --lam 3 --d $d --dps 90 --ygl 24
done
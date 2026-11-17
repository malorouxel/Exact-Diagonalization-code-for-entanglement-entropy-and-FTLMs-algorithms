#!/bin/bash
# Check the serial and MPI builds on small chains.
# Run it from the repository root with make test.
# Set MPIRUN to change the MPI launcher. The default is mpirun -np.

MPIRUN=${MPIRUN:-"mpirun -np"}
ROOT=$(pwd)
OUT=$ROOT/tests/output
mkdir -p "$OUT"
export OMP_NUM_THREADS=1
FAILED=0

# Compare two numbers with a tolerance.
check(){
	local name=$1 a=$2 b=$3 tol=$4
	if awk -v a="$a" -v b="$b" -v t="$tol" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d<t)}'; then
		echo "PASS  $name"
	else
		echo "FAIL  $name : $a vs $b"
		FAILED=1
	fi
}

# Read the ground state energy of a run.
energy(){ grep "Ground state" "$1/Results.txt" | awk '{print $NF}'; }
entropy(){ grep "Entropy" "$1/Results.txt" | awk '{print $NF}'; }

run_example(){
	local example=$1 folder=$2
	cd "$OUT"
	rm -rf "$folder" "${folder}_serial"
	"$ROOT/bin/lanczos" "$ROOT/examples/$example/InputsParameters.txt" "$ROOT/examples/$example/Interactions.txt" > /dev/null
	mv "$folder" "${folder}_serial"
	for np in 1 2 3; do
		rm -rf "$folder"
		$MPIRUN $np "$ROOT/bin/lanczos_mpi" "$ROOT/examples/$example/InputsParameters.txt" "$ROOT/examples/$example/Interactions.txt" > /dev/null
		mv "$folder" "${folder}_mpi$np"
	done
	cd "$ROOT"
}

echo "AKLT spin 1 chain N = 8 with periodic boundaries"
run_example aklt_spin_one AKLT_N8
check "AKLT exact energy -16/3" "$(energy $OUT/AKLT_N8_serial)" -5.333333333333333 1e-5
for np in 1 2 3; do
	check "AKLT serial and MPI $np ranks give the same energy" "$(energy $OUT/AKLT_N8_serial)" "$(energy $OUT/AKLT_N8_mpi$np)" 1e-10
	check "AKLT serial and MPI $np ranks give the same entropy" "$(entropy $OUT/AKLT_N8_serial)" "$(entropy $OUT/AKLT_N8_mpi$np)" 1e-10
done

echo "Heisenberg spin 1/2 chain N = 12 with periodic boundaries"
run_example heisenberg_spin_half Heisenberg_N12
check "Heisenberg reference energy" "$(energy $OUT/Heisenberg_N12_serial)" -5.387390917 1e-5
for np in 1 2 3; do
	check "Heisenberg serial and MPI $np ranks give the same energy" "$(energy $OUT/Heisenberg_N12_serial)" "$(energy $OUT/Heisenberg_N12_mpi$np)" 1e-10
done

echo "Entropy of a saved ground state"
GS=$OUT/Heisenberg_N12_serial/GroundState.bin
S_SERIAL=$("$ROOT/bin/entropy" "$GS" 0,2,4,6,8,10 2> /dev/null | awk '{print $NF}')
S_MPI=$($MPIRUN 2 "$ROOT/bin/entropy_mpi" "$GS" 0,2,4,6,8,10 2> /dev/null | awk '{print $NF}')
S_HALF=$("$ROOT/bin/entropy" "$GS" 0-5 2> /dev/null | awk '{print $NF}')
check "Entropy of sites 0-5 matches the main program" "$S_HALF" "$(entropy $OUT/Heisenberg_N12_serial)" 1e-10
check "Entropy of a non contiguous subsystem serial and MPI" "$S_SERIAL" "$S_MPI" 1e-10
S_ONE=$("$ROOT/bin/entropy" "$GS" 3 2> /dev/null | awk '{print $NF}')
check "One site of a singlet ground state has entropy ln 2" "$S_ONE" 0.6931471805599453 1e-6

if [ $FAILED -eq 0 ]; then
	echo "All tests passed."
else
	echo "Some tests failed."
	exit 1
fi

#include "indexing.h"
#include "tables.h"

// Return the spin 1/2 state with a given index in the sector of fixed charge.
// The charge is the number of up spins.
// States are ranked with the combinatorial number system.
unsigned long long SpinHalfIndexToState(unsigned long long index, int Charge, int N){
	unsigned long long state = 0;
	unsigned long long remainder = index;
	int charge_left = Charge;
	for (int i = 0; i < N; i++){
		unsigned long long count = Binomial[charge_left][N - i - 1];
		if (remainder >= count){
			remainder -= count;
			charge_left--;
			state += Pow2[N - i - 1];
		}
	}
	return state;
}

// Return the index of a spin 1/2 state in its sector of fixed charge.
// Only the up spins contribute. We jump from one up spin to the next.
unsigned long long SpinHalfStateToIndex(unsigned long long state, int N){
	(void)N;
	unsigned long long index = 0;
	int charge = 0;
	while (state){
		int i = __builtin_ctzll(state);
		charge++;
		index += Binomial[charge][i];
		state &= state - 1;
	}
	return index;
}

// Return the spin 1 state with a given index in the sector of fixed charge.
// A spin 1 state is split in two spin 1/2 states.
// The first one marks the sites in state 2.
// The second one gives the state 0 or 1 of the other sites.
// States are sorted by number of sites in state 2 first.
unsigned long long SpinOneIndexToState(unsigned long long index, int Charge, int N){
	long long remainder = (long long)index;
	unsigned long long index_in_block = 0;
	int number_plus = 0;
	while (remainder >= 0){
		index_in_block = (unsigned long long)remainder;
		remainder -= (long long)(Binomial[number_plus][N] * Binomial[Charge - 2 * number_plus][N - number_plus]);
		number_plus++;
	}
	number_plus--;

	unsigned long long block_size = Binomial[Charge - 2 * number_plus][N - number_plus];
	unsigned long long state_plus = SpinHalfIndexToState(index_in_block / block_size, number_plus, N);
	unsigned long long state_rest = SpinHalfIndexToState(index_in_block % block_size, Charge - 2 * number_plus, N - number_plus);

	unsigned long long state = 0;
	int j = 0;
	for (int i = 0; i < N; i++){
		if (state_plus & Pow2[i]){
			state += 2 * Pow3[i];
		}
		else{
			if (state_rest & Pow2[j])
				state += Pow3[i];
			j++;
		}
	}
	return state;
}

// Return the index of a spin 1 state in its sector of fixed charge.
// This is the inverse of SpinOneIndexToState.
unsigned long long SpinOneStateToIndex(unsigned long long state, int Charge, int N){
	unsigned long long quotient = state;
	unsigned long long state_plus = 0;
	unsigned long long state_rest = 0;
	int number_plus = 0;
	int j = 0;
	for (int i = 0; i < N; i++){
		int digit = quotient % 3;
		quotient /= 3;
		if (digit == 2){
			state_plus += Pow2[i];
			number_plus++;
		}
		else{
			state_rest += digit * Pow2[j];
			j++;
		}
	}
	unsigned long long index = Binomial[Charge - 2 * number_plus][N - number_plus] * SpinHalfStateToIndex(state_plus, N);
	index += SpinHalfStateToIndex(state_rest, j);
	for (int k = 0; k < number_plus; k++)
		index += Binomial[k][N] * Binomial[Charge - 2 * k][N - k];
	return index;
}

// Return the index of a basis state in its sector.
// Without the Sz symmetry the index is the state itself.
unsigned long long StateToIndex(unsigned long long state, int N, int Sz, int SzSym, int SpinOne){
	if (SzSym == 0)
		return state;
	if (SpinOne == 1)
		return SpinOneStateToIndex(state, Sz, N);
	return SpinHalfStateToIndex(state, N);
}

// Return the basis state with a given index in its sector.
unsigned long long IndexToState(unsigned long long index, int N, int Sz, int SzSym, int SpinOne){
	if (SzSym == 0)
		return index;
	if (SpinOne == 1)
		return SpinOneIndexToState(index, Sz, N);
	return SpinHalfIndexToState(index, Sz, N);
}

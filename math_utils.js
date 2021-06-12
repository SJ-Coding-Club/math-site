// suite of functions that solve math problems we will add to the site
/**
 * I know younger students haven't gotten to stuff like matrix 
 * elimination, derivatives, integrals, etc., but just contribute
 * whatever you feel comfortable, no matter how simple you may think it is,
 * as it all helps.
 * 
 * Functions to add:
 * - Convert Decimal to simplified Fraction
 * - Summation function
 * - System of equations solver (Matrix fun)
 * - Definite integrals
 * - Derivative at a point
 * 
 * Currently Added:
 * - Greatest common denominator
 * - Simplify fractions
 * - least common multiple
 * - Fast sum of 1 to N integers 
 * - Primality test (non-optimal but definitely can be improved)
 * - Find factors of a number (not optimal either, but works decently for most numbers) 
 * - Helper binary search function to add numbers to an array in order as the array 
 *   grows (might remove this, read below).
 * - factorial of n
 */


// Returns greatest common factor of integers a and b
function gcd(a, b) {
	if (a % b == 0)
		return b;
	var remainder = a % b;
	return gcd(b, remainder);
}

// reduces fraction to simplest form
function simplify_fraction(numerator, denominator) {
	var greatest_common_factor = gcd(numerator, denominator);
	return [numerator / greatest_common_factor, denominator / greatest_common_factor]
}

// returns least common multiple of integers a and b
// proof: https://proofwiki.org/wiki/Product_of_GCD_and_LCM
function lcm(a, b) {
	return a * b / gcd(a, b);
}


// Returns sum of all integers from 1 to n
//
// proof:
// 
// 1 +  2    +  3    + ... + n = S
// n + (n-1) + (n-2) + ... + 1 = S (same S as the top, but in backwards order)
// + 
// ----------------------------------
// (1 + n) + (2 + (n-1)) + (3 + (n-2)) + ... + (n+1) = 2S
//
// so we're really just adding (1 + n), n times
// so 2S = (n + 1) * n
//     S = (n + 1) * n / 2
function natural_sum(n) {
	return (n + 1) * n / 2;
}

// non-optimal primality test
// we only check from [2, sqrt(n)] because that's the smallest range that
// we know will include exactly half the pairs of factors of n that exist 
// and we only need to find the 1st member of a pair to prove its not pirme
function is_prime(n) {
	if (n < 2)
		return false;
	for (var i = 2; i <= Math.sqrt(n); i++) {
		if (n % i == 0)
			return false;
	}
	return true;
}

// binary search func to find index to insert numbers to keep array sorted as it grows.
// But do we even need this? it is doing binary search n times on the factor array,
// so it is O(N * logN), while using array.sort() at the end is also O(N * logN)). 
// Guess I'll have to test it.
function index_to_insert(n, array) {
	var lo = 0;
	var hi = array.length - 1;
	while (lo <= hi) {
		var mid = Math.floor((lo + hi) / 2);
		if (array[mid] < n) 
			lo = mid + 1;
		else if (array[mid] == n) // this will never happen in factor finding function
			return mid + 1;
		else 
			hi = mid - 1;
	}
	return lo;
}

// returns all the factors of n, sorted
// might change function to sort at the end rather than inserting each in order
function get_factors(n) {
	var factors = [1, n];
	for (var i = 2; i <= Math.sqrt(n); i++) 
		if (n % i == 0) {
			var insertion_index = index_to_insert(i, factors);
			if (i * i == n) {
				factors.splice(insertion_index, 0, i);
			}
			else {
				factors.splice(insertion_index, 0, i)
				factors.splice(index_to_insert(n / i, factors), 0, n / i);
			}
		}
	return factors;
}

// returns n! (n factorial)
function factorial(n) {
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}

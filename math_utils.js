/**
 * This is a suite of functions that solve math problems we will add to the site.
 * Right now, it's just a bunch of solvers and calculators, but we can also add practice problems 
 * to the site, where random numbers are generated with a hidden solution using these functions,
 * and the user has to solve it and input the right solution.
 * 
 * I know younger students haven't gotten to stuff like matrix 
 * elimination, derivatives, integrals, etc., but just contribute
 * whatever you feel comfortable, no matter how simple you may think it is,
 * as it all helps.
 * 
 * 
 * Functions to add:
 * - Convert Decimal to simplified Fraction
 * - Summation function
 * - Definite integrals
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
 * - Derivative of f(x) at x
 * - matrix functions
 * - solve system of equations with 2 and 3 variables
 * 
 * Additionally, some of these functions have been implemented abstractly 
 * but not yet added to the site.
 */


// Returns greatest common factor of integers a and b
function gcd(a, b) {
	if (b == 0 || a % b == 0)
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


// function root_finder(a, b, c) { 
// 	var d = Math.sqrt(b*b - 4 * a * c);
// 	if (d < 0)
// 		return "No real roots"
// 	return [(-b + d) / (2 * a), (-b - d) / (2 * a)]
// }

// function root_formatter(a, b, c) {
// 	var roots = root_finder(a,b,c);
// 	if (Arrays.isArray((roots)))
// 		do_thing
// 	return roots;
// }



// calculate d/dx of f at x using the definition of the derivative as a limit. 
function derivativeAtX(f, x) {
	var h = 0.0000001;
	return (f(x + h) - f(x)) / h;
}




// rounds n to places decimal places, default = 2
function round(n, places=2) {
	return Math.round(Math.pow(10, places) * n) / Math.pow(10, places);
}






// matrix stuff, should be separated into its own class
// all of its a bit of a mess right now, can definitely improve it tho



// general function for adding/subtracting matrices to reduce code written
function matrix_add_subtract(a, b, isAdding) {
	if (a.length != b.length || a[0].length != b[0].length)
		throw "Matrices must be the same dimensions"
	var c = [];
	for (var i = 0; i < a.length; i++) {
		c.push([]);
		for (var j = 0; j < a[0].length; j++) {
			if (!isAdding)
				c[i][j] = a[i][j] - b[i][j];
			else
				c[i][j] = a[i][j] + b[i][j];
		}
	}
	return c;	
}

// adds two matrices together
function add_matrices(a, b) {
	return matrix_add_subtract(a, b, true);
}

// subtract one matrix from another
function subtract_matrices(a, b) {
	return matrix_add_subtract(a, b, false);
}

// checks that one matrix equals another
function matrices_are_identical(a, b) {
	if (!Array.isArray(a) || !Array.isArray(b) || a.length != b.length)
		return "A and B must be two m x n matrices"
	for (var i = 0; i < a.length; i++) 
		for (var j = 0; j < a[i].length; j++)
			if (a[i][j] != b[i][j])
				return false;
	return true;
}

// multiplying matrices algorithm, yeah i know it's O(n^3) but we can improve it later
function multiply_matrices(a, b) {
	if (!Array.isArray(a) || !Array.isArray(b))
		throw "A and B must be m x n and n x p matrices"
	if (a[0].length != b.length)
		throw "Number of columns of matrix A must be equal to number of rows of matrix B"
	var c = [];
	// if A is a m x n matrix, B will be n x p matrix. Product C will be m x p matrix
	for (var i = 0; i < a.length; i++) {
		c.push([]);
		for (var j = 0; j < b[0].length; j++) {
		  var dot_product = 0;
		  for (var k = 0; k < b.length; k++) {
		    dot_product += a[i][k] * b[k][j];
		  }
			c[i].push(dot_product);
		}
	}
	return c;
}

function multiply_matrix_by_scalar(a, c) {
	for (var i = 0; i < a.length; i++) {
		for (var j = 0; j < a[i].length; j++) {
			a[i][j] *= c;
		}
	}
	return a;
}

// swaps two rows
function swap_rows(matrix, r1_index, r2_index) {
	var t = matrix[r1_index];
	matrix[r1_index] = matrix[r2_index];
	matrix[r2_index] = t;
	return matrix;
}

// generates identity matrix
function identity_matrix(n) {
	var matrix = [];
	for (var i = 0; i < n; i++) {
		matrix.push(new Array(n).fill(0));
		matrix[i][i] = 1;
	}
	return matrix;
}

// hard coded function to find minors for 3x3
function matrix_of_minors_for_3x3(m) {
	if (m.length != 3 || m[0].length != 3)
		throw "Matrix must be 3x3"
	return [
		[m[1][1]*m[2][2] - m[1][2]*m[2][1], m[1][0]*m[2][2] - m[1][2]*m[2][0], m[1][0]*m[2][1] - m[1][1]*m[2][0]],
		[m[0][1]*m[2][2] - m[0][2]*m[2][1], m[0][0]*m[2][2] - m[2][0]*m[0][2], m[0][0]*m[2][1] - m[2][0]*m[0][1]],
		[m[0][1]*m[1][2] - m[0][2]*m[1][1], m[0][0]*m[1][2] - m[0][2]*m[1][0], m[0][0]*m[1][1] - m[0][1]*m[1][0]]
	];
}

function cofactors(matrix) {
	for (var i = 0; i < matrix.length; i++) {
		for (var j = 0; j < matrix[i].length; j++) {
			matrix[i][j] *= Math.pow(-1, i + j);
		}
	}
	return matrix;
}

// general determinant calculation using Laplace expansion
function determinant(matrix) {
	if (matrix.length == 1) {
		return matrix[0][0];
	}
	var det = 0;
	for (var i = 0; i < matrix[0].length; i++) {
	  var subIndex = 0;
		var subMatrix = [];
		for (var j = 1; j < matrix.length; j++) {
			subMatrix.push([])
			for (var k = 0; k < matrix[j].length; k++) {
				if (k == i) {
					continue;
				}
				subMatrix[subIndex].push(matrix[j][k]);
			}
			subIndex++;
		}
		var factor = i % 2 == 0 ? 1 : -1;
		det += factor * matrix[0][i] * determinant(subMatrix);
	}
	return det;
}

function transpose(matrix) {
  var rows = matrix.length;
  var cols = matrix[0].length;
  new_matrix = []
  for (var i = 0; i < cols; i++) {
    new_matrix.push([]);
    for (var j = 0; j < rows; j++) {
      new_matrix[i].push(matrix[j][i]);
    }
  }
  return new_matrix;
}


// hardcoded method for finding inverse of matrix 
function matrix_inverse(matrix) {
	var len = matrix.length;
	if (len != matrix[0].length)
		throw "Must be a square matrix"
	if (len == 2) {
		var a = matrix[0][0];
		var b = matrix[0][1];
		var c = matrix[1][0];
		var d = matrix[1][1];
		var det = a * d - b * c;
		if (det == 0)
			throw "Inverse matrix does not exist: zero determinate"
		return multiply_matrix_by_scalar([[d, -b], [-c, a]], 1 / det);
	} 
	else if (len == 3) { // using minors, cofactors, and adjugate
		var minors = matrix_of_minors_for_3x3(matrix);
		minors = cofactors(minors);
		minors = transpose(minors);
		var det = determinant(matrix);
		if (det == 0)
			throw "Inverse matrix does not exist: zero determinate";
		minors = multiply_matrix_by_scalar(minors, 1 / det);
		return minors
	}

	else {
		id_matrix = identity_matrix(len);
		// implement better strategy
		// https://www.mathsisfun.com/algebra/matrix-inverse-row-operations-gauss-jordan.html
	}

}

// matrix received in form of the coefficients of system
// from: https://www.mathsisfun.com/algebra/systems-linear-equations-matrices.html
// Example:
//  x +  y +  z =  6
//      2y + 5z = -4
// 2x + 5y -  z = 27
// becomes:
// [ 1  1  1 ] [x]   [  6 ]
// [ 0  2  5 ] [y] = [ -4 ]
// [ 2  5 -1 ] [z]   [ 27 ]
// Can be written as AX = B
// Then X = (A^-1)B

 
function system_solver(coef_matrix, right_matrix) {
	var inv_matrix = matrix_inverse(coef_matrix);
	var solution = multiply_matrices(inv_matrix, right_matrix);
	for (var i = 0; i < solution.length; i++) {
		solution[i][0] = round(solution[i], 8);
	}
	return solution; // clean this up
}


function checkNoVSInfinitelyMany(x1, s1, x2, s2) {
	if (gcd(x1, x2) == gcd(s1, s2))
		return "Infinitely Many Solutions";
	else
		return "No Solution";
}

/*

SLOPE INTERCEPT FORM GENERATOR
Returns a slope intercept form equation of a line from 2 coordinates

*/

function slopeInterceptFormGenerator(a, b) {
	let ax = a[0];
	let ay = a[1];
	let bx = b[0];
	let by = b[1];
	let rise = by - ay;
	let run = bx - ax;
	let slope = simplify_fraction(rise, run);
	rise = slope[0];
	run = slope[1];
	let slopeFraction = `${rise}/${run}`; // Don't want to use decimals
	let slopeDecimal = rise / run; // To use to get the y intercept
	if (slopeDecimal % 1 == 0) slopeFraction = slopeDecimal; // To stop things like "1/1" or "6/3" when they could be whole numbers
	if (slopeFraction == 1) slopeFraction = ""; // To stop "1x"
	let yIntercept = ay - (slopeDecimal * ax);
	if (yIntercept == 0) yIntercept = "";
	else yIntercept = ` + ${yIntercept}`; // To stop "+ 0"
	let final = `y = ${slopeFraction}x${yIntercept}`;
	if (slopeDecimal == 0) {
		if (yIntercept.includes("+")) final = `y = ${yIntercept.slice(3)}`;
		else final = `y = 0`;
	}
	return final;
}


// scripts to interact with html / DOM stuff


function val(id) {
	return parseFloat(document.getElementById(id).value);
}

function parse_gcd() {
	document.getElementById('gcd_answer').innerHTML = gcd(val("gcd_input_a"), val("gcd_input_b"));
}

function parse_lcm() {
	document.getElementById('lcm_answer').innerHTML = lcm(val("lcm_input_a"), val("lcm_input_b"));
}

function parse_factors() {
	document.getElementById('factors_answer').innerHTML = get_factors(val("factors_input")).join(", ");
}

function parse_natural_sum() {
	document.getElementById('natural_sum_answer').innerHTML = natural_sum(val("natural_sum_input"));	
}

function parse_primality_test() {
	var primality = is_prime(val("primality_test_input"));
	document.getElementById('primality_test_answer').innerHTML = primality ? "Prime" : "Not Prime"; 
}


function parse_simplify_fraction() {
	var simplified = simplify_fraction(val("numerator"), val("denominator"));
	document.getElementById('simplified_numerator').innerHTML = simplified[0] + " / ";
	document.getElementById('simplified_denominator').innerHTML = simplified[1];
}


function parse_factorial() {
	document.getElementById('factorial_answer').innerHTML = factorial(val("factorial_input"));
}


function parse_solver_inputs(isTrivariate) {

	if (isTrivariate) {
		var coefficient_matrix = [
			[val("tx1"), val("ty1"), val("tz1")],
			[val("tx2"), val("ty2"), val("tz2")],
			[val("tx3"), val("ty3"), val("tz3")]];
		var right_matrix = [
			[val("ts1")], [val("ts2")], [val("ts3")]
		];
		var solution_message;
		try {
			var solution = system_solver(coefficient_matrix, right_matrix);
			solution_message = "x = " + solution[0][0] + "<br> y = " + solution[1][0] + "<br> z = " + solution[2][0];
		}
		catch(zeroDeterminantError) {
			solution_message = checkNoVSInfinitelyMany(val("tx1"), val("ts1"), val("tx2"), val("ts2"))
		}
		document.getElementById('trivariate_answer').innerHTML = solution_message;

	}
	else {
		var coefficient_matrix = [
			[val("bx1"), val("by1")],
			[val("bx2"), val("by2")]
		];
		var right_matrix = [
			[val("bs1")], [val("bs2")]
		];
		var solution_message;
		try {
			var solution = system_solver(coefficient_matrix, right_matrix);
			solution_message = "x = " + solution[0][0] + "<br> y = " + solution[1][0];
		}
		catch(zeroDeterminantError) {
			solution_message = checkNoVSInfinitelyMany(val("bx1"), val("bs1"), val("bx2"), val("bs2"))
		}
		document.getElementById('bivariate_answer').innerHTML = solution_message;

	}
}

// shows selected div, hides unselected ones
function showDiv(selection) {
	document.getElementById(selection.value).style.display = 'block';
	var functions = document.getElementById('math_functions');
	for (i = 1; i < functions.length; i++) {
		if (functions[i].value != selection.value) {
			document.getElementById(functions[i].value).style.display = 'none';
		} 
	}
}

function parse_slope_intercept() {
	document.getElementById('slope_intercept_answer').innerHTML = slopeInterceptFormGenerator([val("ax"), val("ay")], [val("bx"), val("by")]);
}

REAL slamch(CMACH) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  String CMACH;
  // ..

// =====================================================================

  // .. Parameters ..
  REAL ONE, ZERO;
  const ONE = 1.0, ZERO = 0.0;
  // ..
  // .. Local Scalars ..
  REAL RND, EPS, SFMIN, SMALL, RMACH;
  // ..
  // .. External Functions ..
  //- bool               lsame;
  // EXTERNAL lsame
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
  // ..
  // .. Executable Statements ..

  // Assume rounding, not chopping. Always.

  RND = ONE;

  if (ONE == RND) {
    EPS = EPSILON(ZERO) * 0.5;
  } else {
    EPS = EPSILON(ZERO);
  }

  if (lsame(CMACH, 'E')) {
    RMACH = EPS;
  } else if (lsame(CMACH, 'S')) {
    SFMIN = TINY(ZERO);
    SMALL = ONE / HUGE(ZERO);
    if (SMALL >= SFMIN) {
      // Use SMALL plus a bit, to avoid the possibility of rounding
      // causing overflow when computing  1/sfmin.

      SFMIN = SMALL * (ONE + EPS);
    }
    RMACH = SFMIN;
  } else if (lsame(CMACH, 'B')) {
    RMACH = RADIX(ZERO);
  } else if (lsame(CMACH, 'P')) {
    RMACH = EPS * RADIX(ZERO);
  } else if (lsame(CMACH, 'N')) {
    RMACH = DIGITS(ZERO);
  } else if (lsame(CMACH, 'R')) {
    RMACH = RND;
  } else if (lsame(CMACH, 'M')) {
    RMACH = MINEXPONENT(ZERO);
  } else if (lsame(CMACH, 'U')) {
    RMACH = tiny(zero);
  } else if (lsame(CMACH, 'L')) {
    RMACH = MAXEXPONENT(ZERO);
  } else if (lsame(CMACH, 'O')) {
    RMACH = HUGE(ZERO);
  } else {
    RMACH = ZERO;
  }

  SLAMCH = RMACH;
  return;
}

// ***********************************************************************
// > \brief \b SLAMC3
// > \details
// > \b Purpose:
// > \verbatim
// > SLAMC3  is intended to force  A  and  B  to be stored prior to doing
// > the addition of  A  and  B ,  for use in situations where optimizers
// > might hold one of these in a register.
// > \endverbatim
// > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
// >
// > \ingroup lamc3
// >
// > \param[in] A
// > \verbatim
// > \endverbatim
// >
// > \param[in] B
// > \verbatim
// >          The values A and B.
// > \endverbatim
// >

REAL slamc3(A, B) {
// -- LAPACK auxiliary routine --
  // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

  // .. Scalar Arguments ..
  REAL A, B;
  // ..
// =====================================================================

  // .. Executable Statements ..

  SLAMC3 = A + B;

  return;
}

// ***********************************************************************

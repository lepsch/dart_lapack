// > \ingroup lamch

// =====================================================================
      double dlamch(CMACH ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             CMACH;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      double             RND, EPS, SFMIN, SMALL, RMACH;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
      // ..
      // .. Executable Statements ..


      // Assume rounding, not chopping. Always.

      RND = ONE;

      if ( ONE == RND ) {
         EPS = EPSILON(ZERO) * 0.5;
      } else {
         EPS = EPSILON(ZERO);
      }

      if ( LSAME( CMACH, 'E' ) ) {
         RMACH = EPS;
      } else if ( LSAME( CMACH, 'S' ) ) {
         SFMIN = TINY(ZERO);
         SMALL = ONE / HUGE(ZERO);
         if ( SMALL >= SFMIN ) {

            // Use SMALL plus a bit, to avoid the possibility of rounding
            // causing overflow when computing  1/sfmin.

            SFMIN = SMALL*( ONE+EPS );
         }
         RMACH = SFMIN;
      } else if ( LSAME( CMACH, 'B' ) ) {
         RMACH = RADIX(ZERO);
      } else if ( LSAME( CMACH, 'P' ) ) {
         RMACH = EPS * RADIX(ZERO);
      } else if ( LSAME( CMACH, 'N' ) ) {
         RMACH = DIGITS(ZERO);
      } else if ( LSAME( CMACH, 'R' ) ) {
         RMACH = RND;
      } else if ( LSAME( CMACH, 'M' ) ) {
         RMACH = MINEXPONENT(ZERO);
      } else if ( LSAME( CMACH, 'U' ) ) {
         RMACH = tiny(zero);
      } else if ( LSAME( CMACH, 'L' ) ) {
         RMACH = MAXEXPONENT(ZERO);
      } else if ( LSAME( CMACH, 'O' ) ) {
         RMACH = HUGE(ZERO);
      } else {
         RMACH = ZERO;
      }

      DLAMCH = RMACH;
      return;
      }

// ***********************************************************************
// > \brief \b DLAMC3
// > \details
// > \b Purpose:
// > \verbatim
// > DLAMC3  is intended to force  A  and  B  to be stored prior to doing
// > the addition of  A  and  B ,  for use in situations where optimizers
// > might hold one of these in a register.
// > \endverbatim
// > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
// > \param[in] A
// > \verbatim
// >          A is a DOUBLE PRECISION
// > \endverbatim
// >
// > \param[in] B
// > \verbatim
// >          B is a DOUBLE PRECISION
// >          The values A and B.
// > \endverbatim
// >
// > \ingroup lamc3
// >
      double dlamc3(A, B ) {

// -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

      // .. Scalar Arguments ..
      double             A, B;
      // ..
// =====================================================================

      // .. Executable Statements ..

      DLAMC3 = A + B;

      return;
      }

// ***********************************************************************

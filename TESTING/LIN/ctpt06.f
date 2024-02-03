      SUBROUTINE CTPT06( RCOND, RCONDC, UPLO, DIAG, N, AP, RWORK, RAT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                N;
      REAL               RAT, RCOND, RCONDC
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      REAL               ANORM, BIGNUM, EPS, RMAX, RMIN
      // ..
      // .. External Functions ..
      REAL               CLANTP, SLAMCH
      // EXTERNAL CLANTP, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' )
      RMAX = MAX( RCOND, RCONDC )
      RMIN = MIN( RCOND, RCONDC )

      // Do the easy cases first.

      if ( RMIN < ZERO ) {

         // Invalid value for RCOND or RCONDC, return 1/EPS.

         RAT = ONE / EPS

      } else if ( RMIN > ZERO ) {

         // Both estimates are positive, return RMAX/RMIN - 1.

         RAT = RMAX / RMIN - ONE

      } else if ( RMAX == ZERO ) {

         // Both estimates zero.

         RAT = ZERO

      } else {

         // One estimate is zero, the other is non-zero.  If the matrix is
         // ill-conditioned, return the nonzero estimate multiplied by
         // 1/EPS; if the matrix is badly scaled, return the nonzero
         // estimate multiplied by BIGNUM/TMAX, where TMAX is the maximum
         // element in absolute value in A.

         BIGNUM = ONE / SLAMCH( 'Safe minimum' )
         ANORM = CLANTP( 'M', UPLO, DIAG, N, AP, RWORK )

         RAT = RMAX*( MIN( BIGNUM / MAX( ONE, ANORM ), ONE / EPS ) )
      }

      RETURN

      // End of CTPT06

      }

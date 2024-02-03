      SUBROUTINE CTRT06( RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK, RAT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                LDA, N;
      REAL               RAT, RCOND, RCONDC
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               ANORM, BIGNUM, EPS, RMAX, RMIN
      // ..
      // .. External Functions ..
      REAL               CLANTR, SLAMCH
      // EXTERNAL CLANTR, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'Epsilon' )
      RMAX = MAX( RCOND, RCONDC )
      RMIN = MIN( RCOND, RCONDC )

      // Do the easy cases first.

      IF( RMIN.LT.ZERO ) THEN

         // Invalid value for RCOND or RCONDC, return 1/EPS.

         RAT = ONE / EPS

      ELSE IF( RMIN.GT.ZERO ) THEN

         // Both estimates are positive, return RMAX/RMIN - 1.

         RAT = RMAX / RMIN - ONE

      ELSE IF( RMAX.EQ.ZERO ) THEN

         // Both estimates zero.

         RAT = ZERO

      ELSE

         // One estimate is zero, the other is non-zero.  If the matrix is
         // ill-conditioned, return the nonzero estimate multiplied by
         // 1/EPS; if the matrix is badly scaled, return the nonzero
         // estimate multiplied by BIGNUM/TMAX, where TMAX is the maximum
         // element in absolute value in A.

         BIGNUM = ONE / SLAMCH( 'Safe minimum' )
         ANORM = CLANTR( 'M', UPLO, DIAG, N, N, A, LDA, RWORK )

         RAT = RMAX*( MIN( BIGNUM / MAX( ONE, ANORM ), ONE / EPS ) )
      END IF

      RETURN

      // End of CTRT06

      }

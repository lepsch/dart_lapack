      SUBROUTINE ZTRT06( RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK, RAT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                LDA, N;
      double             RAT, RCOND, RCONDC;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             ANORM, BIGNUM, EPS, RMAX, RMIN;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANTR;
      // EXTERNAL DLAMCH, ZLANTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )
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

         BIGNUM = ONE / DLAMCH( 'Safe minimum' )
         ANORM = ZLANTR( 'M', UPLO, DIAG, N, N, A, LDA, RWORK )

         RAT = RMAX*( MIN( BIGNUM / MAX( ONE, ANORM ), ONE / EPS ) )
      }

      RETURN

      // End of ZTRT06

      }

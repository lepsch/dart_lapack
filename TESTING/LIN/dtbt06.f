      SUBROUTINE DTBT06( RCOND, RCONDC, UPLO, DIAG, N, KD, AB, LDAB, WORK, RAT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                KD, LDAB, N;
      double             RAT, RCOND, RCONDC;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             ANORM, BIGNUM, EPS, RMAX, RMIN, SMLNUM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANTB;
      // EXTERNAL DLAMCH, DLANTB
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )
      RMAX = MAX( RCOND, RCONDC )
      RMIN = MIN( RCOND, RCONDC )

      // Do the easy cases first.

      if ( RMIN.LT.ZERO ) {

         // Invalid value for RCOND or RCONDC, return 1/EPS.

         RAT = ONE / EPS

      } else if ( RMIN.GT.ZERO ) {

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

         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         ANORM = DLANTB( 'M', UPLO, DIAG, N, KD, AB, LDAB, WORK )

         RAT = RMAX*( MIN( BIGNUM / MAX( ONE, ANORM ), ONE / EPS ) )
      }

      RETURN

      // End of DTBT06

      }

      void ctbt06(final int RCOND, final int RCONDC, final int UPLO, final int DIAG, final int N, final int KD, final Matrix<double> AB, final int LDAB, final Array<double> RWORK, final int RAT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                KD, LDAB, N;
      double               RAT, RCOND, RCONDC;
      double               RWORK( * );
      Complex            AB( LDAB, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ANORM, BIGNUM, EPS, RMAX, RMIN;
      // ..
      // .. External Functions ..
      //- REAL               CLANTB, SLAMCH;
      // EXTERNAL CLANTB, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      EPS = SLAMCH( 'Epsilon' );
      RMAX = max( RCOND, RCONDC );
      RMIN = min( RCOND, RCONDC );

      // Do the easy cases first.

      if ( RMIN < ZERO ) {

         // Invalid value for RCOND or RCONDC, return 1/EPS.

         RAT = ONE / EPS;

      } else if ( RMIN > ZERO ) {

         // Both estimates are positive, return RMAX/RMIN - 1.

         RAT = RMAX / RMIN - ONE;

      } else if ( RMAX == ZERO ) {

         // Both estimates zero.

         RAT = ZERO;

      } else {

         // One estimate is zero, the other is non-zero.  If the matrix is
         // ill-conditioned, return the nonzero estimate multiplied by
         // 1/EPS; if the matrix is badly scaled, return the nonzero
         // estimate multiplied by BIGNUM/TMAX, where TMAX is the maximum
         // element in absolute value in A.

         BIGNUM = ONE / SLAMCH( 'Safe minimum' );
         ANORM = CLANTB( 'M', UPLO, DIAG, N, KD, AB, LDAB, RWORK );

         RAT = RMAX*( min( BIGNUM / max( ONE, ANORM ), ONE / EPS ) );
      }

      }

      void dtrt06(final int RCOND, final int RCONDC, final int UPLO, final int DIAG, final int N, final Matrix<double> A_, final int LDA, final Array<double> WORK_, final int RAT,) {
  final A = A_.having();
  final WORK = WORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                LDA, N;
      double             RAT, RCOND, RCONDC;
      double             A( LDA, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ANORM, BIGNUM, EPS, RMAX, RMIN, SMLNUM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANTR;
      // EXTERNAL DLAMCH, DLANTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      EPS = dlamch( 'Epsilon' );
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

         SMLNUM = dlamch( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
         ANORM = DLANTR( 'M', UPLO, DIAG, N, N, A, LDA, WORK );

         RAT = RMAX*( min( BIGNUM / max( ONE, ANORM ), ONE / EPS ) );
      }

      }

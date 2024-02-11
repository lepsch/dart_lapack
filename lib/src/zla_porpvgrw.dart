      double zla_porpvgrw(final int UPLO, final int NCOLS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final Array<double> WORK,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                NCOLS, LDA, LDAF;
      Complex         A( LDA, * ), AF( LDAF, * );
      double             WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             AMAX, UMAX, RPVGRW;
      bool               UPPER;
      Complex         ZDUM;
      // ..
      // .. External Functions ..
      // EXTERNAL lsame
      bool               lsame;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      UPPER = lsame( 'Upper', UPLO );

      // DPOTRF will have factored only the NCOLSxNCOLS leading submatrix,
      // so we restrict the growth search to that submatrix and use only
      // the first 2*NCOLS workspace entries.

      RPVGRW = 1.0;
      for (I = 1; I <= 2*NCOLS; I++) {
         WORK[I] = 0.0;
      }

      // Find the max magnitude entry of each column.

      if ( UPPER ) {
         for (J = 1; J <= NCOLS; J++) {
            for (I = 1; I <= J; I++) {
               WORK[NCOLS+J] = max( CABS1( A( I, J ) ), WORK( NCOLS+J ) );
            }
         }
      } else {
         for (J = 1; J <= NCOLS; J++) {
            for (I = J; I <= NCOLS; I++) {
               WORK[NCOLS+J] = max( CABS1( A( I, J ) ), WORK( NCOLS+J ) );
            }
         }
      }

      // Now find the max magnitude entry of each column of the factor in
      // AF.  No pivoting, so no permutations.

      if ( lsame( 'Upper', UPLO ) ) {
         for (J = 1; J <= NCOLS; J++) {
            for (I = 1; I <= J; I++) {
               WORK[J] = max( CABS1( AF( I, J ) ), WORK( J ) );
            }
         }
      } else {
         for (J = 1; J <= NCOLS; J++) {
            for (I = J; I <= NCOLS; I++) {
               WORK[J] = max( CABS1( AF( I, J ) ), WORK( J ) );
            }
         }
      }

      // Compute the *inverse* of the max element growth factor.  Dividing
      // by zero would imply the largest entry of the factor's column is
      // zero.  Than can happen when either the column of A is zero or
      // massive pivots made the factor underflow to zero.  Neither counts
      // as growth in itself, so simply ignore terms with zero
      // denominators.

      if ( lsame( 'Upper', UPLO ) ) {
         for (I = 1; I <= NCOLS; I++) {
            UMAX = WORK( I );
            AMAX = WORK( NCOLS+I );
            if ( UMAX /= 0.0 ) {
               RPVGRW = min( AMAX / UMAX, RPVGRW );
            }
         }
      } else {
         for (I = 1; I <= NCOLS; I++) {
            UMAX = WORK( I );
            AMAX = WORK( NCOLS+I );
            if ( UMAX /= 0.0 ) {
               RPVGRW = min( AMAX / UMAX, RPVGRW );
            }
         }
      }

      ZLA_PORPVGRW = RPVGRW;
      }

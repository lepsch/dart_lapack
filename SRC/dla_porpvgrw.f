      double           FUNCTION DLA_PORPVGRW( UPLO, NCOLS, A, LDA, AF, LDAF, WORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                NCOLS, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AF( LDAF, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             AMAX, UMAX, RPVGRW;
      bool               UPPER;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME
      bool               LSAME;
      // ..
      // .. Executable Statements ..

      UPPER = LSAME( 'Upper', UPLO )

      // DPOTRF will have factored only the NCOLSxNCOLS leading submatrix,
      // so we restrict the growth search to that submatrix and use only
      // the first 2*NCOLS workspace entries.

      RPVGRW = 1.0D+0
      DO I = 1, 2*NCOLS
         WORK( I ) = 0.0D+0
      END DO

      // Find the max magnitude entry of each column.

      if ( UPPER ) {
         DO J = 1, NCOLS
            DO I = 1, J
               WORK( NCOLS+J ) = MAX( ABS( A( I, J ) ), WORK( NCOLS+J ) )
            END DO
         END DO
      } else {
         DO J = 1, NCOLS
            DO I = J, NCOLS
               WORK( NCOLS+J ) = MAX( ABS( A( I, J ) ), WORK( NCOLS+J ) )
            END DO
         END DO
      }

      // Now find the max magnitude entry of each column of the factor in
      // AF.  No pivoting, so no permutations.

      if ( LSAME( 'Upper', UPLO ) ) {
         DO J = 1, NCOLS
            DO I = 1, J
               WORK( J ) = MAX( ABS( AF( I, J ) ), WORK( J ) )
            END DO
         END DO
      } else {
         DO J = 1, NCOLS
            DO I = J, NCOLS
               WORK( J ) = MAX( ABS( AF( I, J ) ), WORK( J ) )
            END DO
         END DO
      }

      // Compute the *inverse* of the max element growth factor.  Dividing
      // by zero would imply the largest entry of the factor's column is
      // zero.  Than can happen when either the column of A is zero or
      // massive pivots made the factor underflow to zero.  Neither counts
      // as growth in itself, so simply ignore terms with zero
      // denominators.

      if ( LSAME( 'Upper', UPLO ) ) {
         DO I = 1, NCOLS
            UMAX = WORK( I )
            AMAX = WORK( NCOLS+I )
            if ( UMAX /= 0.0D+0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         END DO
      } else {
         DO I = 1, NCOLS
            UMAX = WORK( I )
            AMAX = WORK( NCOLS+I )
            if ( UMAX /= 0.0D+0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         END DO
      }

      DLA_PORPVGRW = RPVGRW

      // End of DLA_PORPVGRW

      }

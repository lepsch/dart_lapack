      REAL FUNCTION SLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, INFO, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), AF( LDAF, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                NCOLS, I, J, K, KP;
      REAL               AMAX, UMAX, RPVGRW, TMP
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
      if ( INFO.EQ.0 ) {
         if ( UPPER ) {
            NCOLS = 1
         } else {
            NCOLS = N
         }
      } else {
         NCOLS = INFO
      }

      RPVGRW = 1.0
      for (I = 1; I <= 2*N; I++) {
         WORK( I ) = 0.0
      END DO

      // Find the max magnitude entry of each column of A.  Compute the max
      // for all N columns so we can apply the pivot permutation while
      // looping below.  Assume a full factorization is the common case.

      if ( UPPER ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               WORK( N+I ) = MAX( ABS( A( I, J ) ), WORK( N+I ) )
               WORK( N+J ) = MAX( ABS( A( I, J ) ), WORK( N+J ) )
            END DO
         END DO
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               WORK( N+I ) = MAX( ABS( A( I, J ) ), WORK( N+I ) )
               WORK( N+J ) = MAX( ABS( A( I, J ) ), WORK( N+J ) )
            END DO
         END DO
      }

      // Now find the max magnitude entry of each column of U or L.  Also
      // permute the magnitudes of A above so they're in the same order as
      // the factor.

      // The iteration orders and permutations were copied from ssytrs.
      // Calls to SSWAP would be severe overkill.

      if ( UPPER ) {
         K = N
         DO WHILE ( K .LT. NCOLS .AND. K.GT.0 )
            if ( IPIV( K ).GT.0 ) {
               // 1x1 pivot
               KP = IPIV( K )
               if ( KP .NE. K ) {
                  TMP = WORK( N+K )
                  WORK( N+K ) = WORK( N+KP )
                  WORK( N+KP ) = TMP
               }
               for (I = 1; I <= K; I++) {
                  WORK( K ) = MAX( ABS( AF( I, K ) ), WORK( K ) )
               END DO
               K = K - 1
            } else {
               // 2x2 pivot
               KP = -IPIV( K )
               TMP = WORK( N+K-1 )
               WORK( N+K-1 ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               for (I = 1; I <= K-1; I++) {
                  WORK( K ) = MAX( ABS( AF( I, K ) ), WORK( K ) )
                  WORK( K-1 ) = MAX( ABS( AF( I, K-1 ) ), WORK( K-1 ) )
               END DO
               WORK( K ) = MAX( ABS( AF( K, K ) ), WORK( K ) )
               K = K - 2
            }
         END DO
         K = NCOLS
         DO WHILE ( K .LE. N )
            if ( IPIV( K ).GT.0 ) {
               KP = IPIV( K )
               if ( KP .NE. K ) {
                  TMP = WORK( N+K )
                  WORK( N+K ) = WORK( N+KP )
                  WORK( N+KP ) = TMP
               }
               K = K + 1
            } else {
               KP = -IPIV( K )
               TMP = WORK( N+K )
               WORK( N+K ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               K = K + 2
            }
         END DO
      } else {
         K = 1
         DO WHILE ( K .LE. NCOLS )
            if ( IPIV( K ).GT.0 ) {
               // 1x1 pivot
               KP = IPIV( K )
               if ( KP .NE. K ) {
                  TMP = WORK( N+K )
                  WORK( N+K ) = WORK( N+KP )
                  WORK( N+KP ) = TMP
               }
               for (I = K; I <= N; I++) {
                  WORK( K ) = MAX( ABS( AF( I, K ) ), WORK( K ) )
               END DO
               K = K + 1
            } else {
               // 2x2 pivot
               KP = -IPIV( K )
               TMP = WORK( N+K+1 )
               WORK( N+K+1 ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               for (I = K+1; I <= N; I++) {
                  WORK( K ) = MAX( ABS( AF( I, K ) ), WORK( K ) )
                  WORK( K+1 ) = MAX( ABS( AF(I, K+1 ) ), WORK( K+1 ) )
               END DO
               WORK( K ) = MAX( ABS( AF( K, K ) ), WORK( K ) )
               K = K + 2
            }
         END DO
         K = NCOLS
         DO WHILE ( K .GE. 1 )
            if ( IPIV( K ).GT.0 ) {
               KP = IPIV( K )
               if ( KP .NE. K ) {
                  TMP = WORK( N+K )
                  WORK( N+K ) = WORK( N+KP )
                  WORK( N+KP ) = TMP
               }
               K = K - 1
            } else {
               KP = -IPIV( K )
               TMP = WORK( N+K )
               WORK( N+K ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               K = K - 2
            }
         END DO
      }

      // Compute the *inverse* of the max element growth factor.  Dividing
      // by zero would imply the largest entry of the factor's column is
      // zero.  Than can happen when either the column of A is zero or
      // massive pivots made the factor underflow to zero.  Neither counts
      // as growth in itself, so simply ignore terms with zero
      // denominators.

      if ( UPPER ) {
         for (I = NCOLS; I <= N; I++) {
            UMAX = WORK( I )
            AMAX = WORK( N+I )
            if ( UMAX /= 0.0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         END DO
      } else {
         for (I = 1; I <= NCOLS; I++) {
            UMAX = WORK( I )
            AMAX = WORK( N+I )
            if ( UMAX /= 0.0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         END DO
      }

      SLA_SYRPVGRW = RPVGRW

      // End of SLA_SYRPVGRW

      }

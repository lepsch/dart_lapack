      REAL FUNCTION CLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, INFO, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), AF( LDAF, * )
      REAL               WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                NCOLS, I, J, K, KP;
      REAL               AMAX, UMAX, RPVGRW, TMP
      bool               UPPER, LSAME;
      COMPLEX            ZDUM
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, AIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL ( ZDUM ) ) + ABS( AIMAG ( ZDUM ) )
      // ..
      // .. Executable Statements ..

      UPPER = LSAME( 'Upper', UPLO )
      if ( INFO == 0 ) {
         if (UPPER) {
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
      }

      // Find the max magnitude entry of each column of A.  Compute the max
      // for all N columns so we can apply the pivot permutation while
      // looping below.  Assume a full factorization is the common case.

      if ( UPPER ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               WORK( N+I ) = MAX( CABS1( A( I,J ) ), WORK( N+I ) )
               WORK( N+J ) = MAX( CABS1( A( I,J ) ), WORK( N+J ) )
            }
         }
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               WORK( N+I ) = MAX( CABS1( A( I, J ) ), WORK( N+I ) )
               WORK( N+J ) = MAX( CABS1( A( I, J ) ), WORK( N+J ) )
            }
         }
      }

      // Now find the max magnitude entry of each column of U or L.  Also
      // permute the magnitudes of A above so they're in the same order as
      // the factor.

      // The iteration orders and permutations were copied from csytrs.
      // Calls to SSWAP would be severe overkill.

      if ( UPPER ) {
         K = N
         DO WHILE ( K < NCOLS && K.GT.0 )
            if ( IPIV( K ).GT.0 ) {
               // 1x1 pivot
               KP = IPIV( K )
               if ( KP != K ) {
                  TMP = WORK( N+K )
                  WORK( N+K ) = WORK( N+KP )
                  WORK( N+KP ) = TMP
               }
               for (I = 1; I <= K; I++) {
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
               }
               K = K - 1
            } else {
               // 2x2 pivot
               KP = -IPIV( K )
               TMP = WORK( N+K-1 )
               WORK( N+K-1 ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               for (I = 1; I <= K-1; I++) {
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
                  WORK( K-1 ) = MAX( CABS1( AF( I, K-1 ) ), WORK( K-1 ) )
               }
               WORK( K ) = MAX( CABS1( AF( K, K ) ), WORK( K ) )
               K = K - 2
            }
         }
         K = NCOLS
         DO WHILE ( K .LE. N )
            if ( IPIV( K ).GT.0 ) {
               KP = IPIV( K )
               if ( KP != K ) {
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
         }
      } else {
         K = 1
         DO WHILE ( K .LE. NCOLS )
            if ( IPIV( K ).GT.0 ) {
               // 1x1 pivot
               KP = IPIV( K )
               if ( KP != K ) {
                  TMP = WORK( N+K )
                  WORK( N+K ) = WORK( N+KP )
                  WORK( N+KP ) = TMP
               }
               for (I = K; I <= N; I++) {
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
               }
               K = K + 1
            } else {
               // 2x2 pivot
               KP = -IPIV( K )
               TMP = WORK( N+K+1 )
               WORK( N+K+1 ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               for (I = K+1; I <= N; I++) {
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
                  WORK( K+1 ) = MAX( CABS1( AF( I, K+1 ) ) , WORK( K+1 ) )
               }
               WORK(K) = MAX( CABS1( AF( K, K ) ), WORK( K ) )
               K = K + 2
            }
         }
         K = NCOLS
         DO WHILE ( K .GE. 1 )
            if ( IPIV( K ).GT.0 ) {
               KP = IPIV( K )
               if ( KP != K ) {
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
         }
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
         }
      } else {
         for (I = 1; I <= NCOLS; I++) {
            UMAX = WORK( I )
            AMAX = WORK( N+I )
            if ( UMAX /= 0.0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         }
      }

      CLA_HERPVGRW = RPVGRW

      // End of CLA_HERPVGRW

      }

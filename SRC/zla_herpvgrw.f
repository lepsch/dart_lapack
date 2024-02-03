      double           FUNCTION ZLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, INFO, LDA, LDAF;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * )
      double             WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                NCOLS, I, J, K, KP;
      double             AMAX, UMAX, RPVGRW, TMP;
      bool               UPPER, LSAME;
      COMPLEX*16         ZDUM
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE ( ZDUM ) ) + ABS( DIMAG ( ZDUM ) )
      // ..
      // .. Executable Statements ..

      UPPER = LSAME( 'Upper', UPLO )
      if ( INFO.EQ.0 ) {
         if (UPPER) {
            NCOLS = 1
         } else {
            NCOLS = N
         }
      } else {
         NCOLS = INFO
      }

      RPVGRW = 1.0D+0
      DO I = 1, 2*N
         WORK( I ) = 0.0D+0
      END DO

      // Find the max magnitude entry of each column of A.  Compute the max
      // for all N columns so we can apply the pivot permutation while
      // looping below.  Assume a full factorization is the common case.

      if ( UPPER ) {
         DO J = 1, N
            DO I = 1, J
               WORK( N+I ) = MAX( CABS1( A( I,J ) ), WORK( N+I ) )
               WORK( N+J ) = MAX( CABS1( A( I,J ) ), WORK( N+J ) )
            END DO
         END DO
      } else {
         DO J = 1, N
            DO I = J, N
               WORK( N+I ) = MAX( CABS1( A( I, J ) ), WORK( N+I ) )
               WORK( N+J ) = MAX( CABS1( A( I, J ) ), WORK( N+J ) )
            END DO
         END DO
      }

      // Now find the max magnitude entry of each column of U or L.  Also
      // permute the magnitudes of A above so they're in the same order as
     t // he factor.

      // The iteration orders and permutations were copied from zsytrs.
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
               DO I = 1, K
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
               END DO
               K = K - 1
            } else {
               // 2x2 pivot
               KP = -IPIV( K )
               TMP = WORK( N+K-1 )
               WORK( N+K-1 ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               DO I = 1, K-1
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
                  WORK( K-1 ) = MAX( CABS1( AF( I, K-1 ) ), WORK( K-1 ) )
               END DO
               WORK( K ) = MAX( CABS1( AF( K, K ) ), WORK( K ) )
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
               DO I = K, N
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
               END DO
               K = K + 1
            } else {
               // 2x2 pivot
               KP = -IPIV( K )
               TMP = WORK( N+K+1 )
               WORK( N+K+1 ) = WORK( N+KP )
               WORK( N+KP ) = TMP
               DO I = K+1, N
                  WORK( K ) = MAX( CABS1( AF( I, K ) ), WORK( K ) )
                  WORK( K+1 ) = MAX( CABS1( AF( I, K+1 ) ) , WORK( K+1 ) )
               END DO
               WORK(K) = MAX( CABS1( AF( K, K ) ), WORK( K ) )
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
            ENDIF
         END DO
      }

      // Compute the *inverse* of the max element growth factor.  Dividing
      // by zero would imply the largest entry of the factor's column is
      // zero.  Than can happen when either the column of A is zero or
      // massive pivots made the factor underflow to zero.  Neither counts
      // as growth in itself, so simply ignore terms with zero
      // denominators.

      if ( UPPER ) {
         DO I = NCOLS, N
            UMAX = WORK( I )
            AMAX = WORK( N+I )
            if ( UMAX /= 0.0D+0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         END DO
      } else {
         DO I = 1, NCOLS
            UMAX = WORK( I )
            AMAX = WORK( N+I )
            if ( UMAX /= 0.0D+0 ) {
               RPVGRW = MIN( AMAX / UMAX, RPVGRW )
            }
         END DO
      }

      ZLA_HERPVGRW = RPVGRW

      // End of ZLA_HERPVGRW

      }

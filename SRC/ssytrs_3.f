      SUBROUTINE SSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), B( LDB, * ), E( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J, K, KP;
      REAL               AK, AKM1, AKM1K, BK, BKM1, DENOM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSWAP, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SSYTRS_3', -INFO );
         RETURN;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      if ( UPPER ) {

         // Begin Upper

         // Solve A*X = B, where A = U*D*U**T.

         // P**T * B

         // Interchange rows K and IPIV(K) of matrix B in the same order
         // that the formation order of IPIV(I) vector for Upper case.

         // (We can do the simple loop over IPIV with decrement -1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = N, 1, -1;
            KP = ABS( IPIV( K ) );
            if ( KP != K ) {
               sswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

         // Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

         strsm('L', 'U', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // Compute D \ B -> B   [ D \ (U \P**T * B) ]

         I = N;
         DO WHILE ( I >= 1 );
            if ( IPIV( I ) > 0 ) {
               sscal(NRHS, ONE / A( I, I ), B( I, 1 ), LDB );
            } else if ( I > 1 ) {
               AKM1K = E( I );
               AKM1 = A( I-1, I-1 ) / AKM1K;
               AK = A( I, I ) / AKM1K;
               DENOM = AKM1*AK - ONE;
               for (J = 1; J <= NRHS; J++) {
                  BKM1 = B( I-1, J ) / AKM1K;
                  BK = B( I, J ) / AKM1K;
                  B( I-1, J ) = ( AK*BKM1-BK ) / DENOM;
                  B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM;
               }
               I = I - 1;
            }
            I = I - 1;
         }

         // Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]

         strsm('L', 'U', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]

         // Interchange rows K and IPIV(K) of matrix B in reverse order
         // from the formation order of IPIV(I) vector for Upper case.

         // (We can do the simple loop over IPIV with increment 1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = 1, N, 1;
            KP = ABS( IPIV( K ) );
            if ( KP != K ) {
               sswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

      } else {

         // Begin Lower

         // Solve A*X = B, where A = L*D*L**T.

         // P**T * B
         // Interchange rows K and IPIV(K) of matrix B in the same order
         // that the formation order of IPIV(I) vector for Lower case.

         // (We can do the simple loop over IPIV with increment 1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = 1, N, 1;
            KP = ABS( IPIV( K ) );
            if ( KP != K ) {
               sswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

         // Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

         strsm('L', 'L', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // Compute D \ B -> B   [ D \ (L \P**T * B) ]

         I = 1;
         DO WHILE ( I <= N );
            if ( IPIV( I ) > 0 ) {
               sscal(NRHS, ONE / A( I, I ), B( I, 1 ), LDB );
            } else if ( I < N ) {
               AKM1K = E( I );
               AKM1 = A( I, I ) / AKM1K;
               AK = A( I+1, I+1 ) / AKM1K;
               DENOM = AKM1*AK - ONE;
               for (J = 1; J <= NRHS; J++) {
                  BKM1 = B( I, J ) / AKM1K;
                  BK = B( I+1, J ) / AKM1K;
                  B( I, J ) = ( AK*BKM1-BK ) / DENOM;
                  B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM;
               }
               I = I + 1;
            }
            I = I + 1;
         }

         // Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]

         strsm('L', 'L', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]

         // Interchange rows K and IPIV(K) of matrix B in reverse order
         // from the formation order of IPIV(I) vector for Lower case.

         // (We can do the simple loop over IPIV with decrement -1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = N, 1, -1;
            KP = ABS( IPIV( K ) );
            if ( KP != K ) {
               sswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

         // END Lower

      }

      RETURN;

      // End of SSYTRS_3

      }

      SUBROUTINE ZSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0,0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J, K, KP;
      COMPLEX*16         AK, AKM1, AKM1K, BK, BKM1, DENOM
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZSCAL, ZSWAP, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('ZSYTRS_3', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0 .OR. NRHS.EQ.0) RETURN;

      if ( UPPER ) {

         // Begin Upper

         // Solve A*X = B, where A = U*D*U**T.

         // P**T * B

         // Interchange rows K and IPIV(K) of matrix B in the same order
         // that the formation order of IPIV(I) vector for Upper case.

         // (We can do the simple loop over IPIV with decrement -1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = N, 1, -1
            KP = ABS( IPIV( K ) )
            if ( KP.NE.K ) {
               zswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

         // Compute (U \P**T * B) -> B    [ (U \P**T * B) ]

         ztrsm('L', 'U', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // Compute D \ B -> B   [ D \ (U \P**T * B) ]

         I = N
         DO WHILE ( I.GE.1 )
            if ( IPIV( I ).GT.0 ) {
               zscal(NRHS, ONE / A( I, I ), B( I, 1 ), LDB );
            } else if ( I.GT.1 ) {
               AKM1K = E( I )
               AKM1 = A( I-1, I-1 ) / AKM1K
               AK = A( I, I ) / AKM1K
               DENOM = AKM1*AK - ONE
               for (J = 1; J <= NRHS; J++) {
                  BKM1 = B( I-1, J ) / AKM1K
                  BK = B( I, J ) / AKM1K
                  B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                  B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
               }
               I = I - 1
            }
            I = I - 1
         }

         // Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]

         ztrsm('L', 'U', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]

         // Interchange rows K and IPIV(K) of matrix B in reverse order
         // from the formation order of IPIV(I) vector for Upper case.

         // (We can do the simple loop over IPIV with increment 1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = 1, N, 1
            KP = ABS( IPIV( K ) )
            if ( KP.NE.K ) {
               zswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
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

         DO K = 1, N, 1
            KP = ABS( IPIV( K ) )
            if ( KP.NE.K ) {
               zswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

         // Compute (L \P**T * B) -> B    [ (L \P**T * B) ]

         ztrsm('L', 'L', 'N', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // Compute D \ B -> B   [ D \ (L \P**T * B) ]

         I = 1
         DO WHILE ( I.LE.N )
            if ( IPIV( I ).GT.0 ) {
               zscal(NRHS, ONE / A( I, I ), B( I, 1 ), LDB );
            } else if ( I.LT.N ) {
               AKM1K = E( I )
               AKM1 = A( I, I ) / AKM1K
               AK = A( I+1, I+1 ) / AKM1K
               DENOM = AKM1*AK - ONE
               for (J = 1; J <= NRHS; J++) {
                  BKM1 = B( I, J ) / AKM1K
                  BK = B( I+1, J ) / AKM1K
                  B( I, J ) = ( AK*BKM1-BK ) / DENOM
                  B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
               }
               I = I + 1
            }
            I = I + 1
         }

         // Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]

         ztrsm('L', 'L', 'T', 'U', N, NRHS, ONE, A, LDA, B, LDB );

         // P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]

         // Interchange rows K and IPIV(K) of matrix B in reverse order
         // from the formation order of IPIV(I) vector for Lower case.

         // (We can do the simple loop over IPIV with decrement -1,
         // since the ABS value of IPIV(I) represents the row index
         // of the interchange with row i in both 1x1 and 2x2 pivot cases)

         DO K = N, 1, -1
            KP = ABS( IPIV( K ) )
            if ( KP.NE.K ) {
               zswap(NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

         // END Lower

      }

      RETURN

      // End of ZSYTRS_3

      }

      SUBROUTINE ZLAVHE( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT;
      int                J, K, KP;
      COMPLEX*16         D11, D12, D21, D22, T1, T2
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZGERU, ZLACGV, ZSCAL, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) && .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.LSAME( DIAG, 'U' ) && .NOT.LSAME( DIAG, 'N' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('ZLAVHE ', -INFO );
         RETURN
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      NOUNIT = LSAME( DIAG, 'N' )
*------------------------------------------

      // Compute  B := A * B  (No transpose)

*------------------------------------------
      if ( LSAME( TRANS, 'N' ) ) {

         // Compute  B := U*B
         // where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))

         if ( LSAME( UPLO, 'U' ) ) {

         // Loop forward applying the transformations.

            K = 1
            } // 10
            if (K.GT.N) GO TO 30;
            if ( IPIV( K ).GT.0 ) {

               // 1 x 1 pivot block

               // Multiply by the diagonal element if forming U * D.

               if (NOUNIT) CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );

               // Multiply by  P(K) * inv(U(K))  if K > 1.

               if ( K.GT.1 ) {

                  // Apply the transformation.

                  zgeru(K-1, NRHS, ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K + 1
            } else {

               // 2 x 2 pivot block

               // Multiply by the diagonal block if forming U * D.

               if ( NOUNIT ) {
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D12 = A( K, K+1 )
                  D21 = DCONJG( D12 )
                  for (J = 1; J <= NRHS; J++) { // 20
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
                  } // 20
               }

               // Multiply by  P(K) * inv(U(K))  if K > 1.

               if ( K.GT.1 ) {

                  // Apply the transformations.

                  zgeru(K-1, NRHS, ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );
                  zgeru(K-1, NRHS, ONE, A( 1, K+1 ), 1, B( K+1, 1 ), LDB, B( 1, 1 ), LDB );

                  // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) )
                  if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K + 2
            }
            GO TO 10
            } // 30

         // Compute  B := L*B
         // where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .

         } else {

            // Loop backward applying the transformations to B.

            K = N
            } // 40
            if (K.LT.1) GO TO 60;

            // Test the pivot index.  If greater than zero, a 1 x 1
            // pivot was used, otherwise a 2 x 2 pivot was used.

            if ( IPIV( K ).GT.0 ) {

               // 1 x 1 pivot block:

               // Multiply by the diagonal element if forming L * D.

               if (NOUNIT) CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K != N ) {
                  KP = IPIV( K )

                  // Apply the transformation.

                  zgeru(N-K, NRHS, ONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K - 1

            } else {

               // 2 x 2 pivot block:

               // Multiply by the diagonal block if forming L * D.

               if ( NOUNIT ) {
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D21 = A( K, K-1 )
                  D12 = DCONJG( D21 )
                  for (J = 1; J <= NRHS; J++) { // 50
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
                  } // 50
               }

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K != N ) {

                  // Apply the transformation.

                  zgeru(N-K, NRHS, ONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );
                  zgeru(N-K, NRHS, ONE, A( K+1, K-1 ), 1, B( K-1, 1 ), LDB, B( K+1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  KP = ABS( IPIV( K ) )
                  if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K - 2
            }
            GO TO 40
            } // 60
         }
*--------------------------------------------------

      // Compute  B := A^H * B  (conjugate transpose)

*--------------------------------------------------
      } else {

         // Form  B := U^H*B
         // where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
         // and   U^H = inv(U^H(1))*P(1)* ... *inv(U^H(m))*P(m)

         if ( LSAME( UPLO, 'U' ) ) {

            // Loop backward applying the transformations.

            K = N
            } // 70
            if (K.LT.1) GO TO 90;

            // 1 x 1 pivot block.

            if ( IPIV( K ).GT.0 ) {
               if ( K.GT.1 ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation
                     // y = y - B' conjg(x),
                  // where x is a column of A and y is a row of B.

                  zlacgv(NRHS, B( K, 1 ), LDB );
                  zgemv('Conjugate', K-1, NRHS, ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB );
                  zlacgv(NRHS, B( K, 1 ), LDB );
               }
               if (NOUNIT) CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );
               K = K - 1

            // 2 x 2 pivot block.

            } else {
               if ( K.GT.2 ) {

                  // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) )
                  if (KP != K-1) CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformations
                     // y = y - B' conjg(x),
                  // where x is a block column of A and y is a block
                  // row of B.

                  zlacgv(NRHS, B( K, 1 ), LDB );
                  zgemv('Conjugate', K-2, NRHS, ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB );
                  zlacgv(NRHS, B( K, 1 ), LDB );

                  zlacgv(NRHS, B( K-1, 1 ), LDB );
                  zgemv('Conjugate', K-2, NRHS, ONE, B, LDB, A( 1, K-1 ), 1, ONE, B( K-1, 1 ), LDB );
                  zlacgv(NRHS, B( K-1, 1 ), LDB );
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D12 = A( K-1, K )
                  D21 = DCONJG( D12 )
                  for (J = 1; J <= NRHS; J++) { // 80
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
                  } // 80
               }
               K = K - 2
            }
            GO TO 70
            } // 90

         // Form  B := L^H*B
         // where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
         // and   L^H = inv(L^H(m))*P(m)* ... *inv(L^H(1))*P(1)

         } else {

            // Loop forward applying the L-transformations.

            K = 1
            } // 100
            if (K.GT.N) GO TO 120;

            // 1 x 1 pivot block

            if ( IPIV( K ).GT.0 ) {
               if ( K.LT.N ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  if (KP != K) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  zlacgv(NRHS, B( K, 1 ), LDB );
                  zgemv('Conjugate', N-K, NRHS, ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB );
                  zlacgv(NRHS, B( K, 1 ), LDB );
               }
               if (NOUNIT) CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );
               K = K + 1

            // 2 x 2 pivot block.

            } else {
               if ( K.LT.N-1 ) {

               // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) )
                  if (KP != K+1) CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  zlacgv(NRHS, B( K+1, 1 ), LDB );
                  zgemv('Conjugate', N-K-1, NRHS, ONE, B( K+2, 1 ), LDB, A( K+2, K+1 ), 1, ONE, B( K+1, 1 ), LDB );
                  zlacgv(NRHS, B( K+1, 1 ), LDB );

                  zlacgv(NRHS, B( K, 1 ), LDB );
                  zgemv('Conjugate', N-K-1, NRHS, ONE, B( K+2, 1 ), LDB, A( K+2, K ), 1, ONE, B( K, 1 ), LDB );
                  zlacgv(NRHS, B( K, 1 ), LDB );
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D21 = A( K+1, K )
                  D12 = DCONJG( D21 )
                  for (J = 1; J <= NRHS; J++) { // 110
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
                  } // 110
               }
               K = K + 2
            }
            GO TO 100
            } // 120
         }

      }
      RETURN

      // End of ZLAVHE

      }

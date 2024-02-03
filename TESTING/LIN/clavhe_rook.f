      SUBROUTINE CLAVHE_ROOK( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT;
      int                J, K, KP;
      COMPLEX            D11, D12, D21, D22, T1, T2
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERU, CLACGV, CSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX
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
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('CLAVHE_ROOK ', -INFO );
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
            if (K > N) GO TO 30;
            if ( IPIV( K ) > 0 ) {

               // 1 x 1 pivot block

               // Multiply by the diagonal element if forming U * D.

               if (NOUNIT) CALL CSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );

               // Multiply by  P(K) * inv(U(K))  if K > 1.

               if ( K > 1 ) {

                  // Apply the transformation.

                  cgeru(K-1, NRHS, CONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K + 1
            } else {

               // 2 x 2 pivot block

               // Multiply by the diagonal block if forming U * D.

               if ( NOUNIT ) {
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D12 = A( K, K+1 )
                  D21 = CONJG( D12 )
                  for (J = 1; J <= NRHS; J++) { // 20
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
                  } // 20
               }

               // Multiply by  P(K) * inv(U(K))  if K > 1.

               if ( K > 1 ) {

                  // Apply the transformations.

                  cgeru(K-1, NRHS, CONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );
                  cgeru(K-1, NRHS, CONE, A( 1, K+1 ), 1, B( K+1, 1 ), LDB, B( 1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  // Swap the first of pair with IMAXth

                  KP = ABS( IPIV( K ) )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // NOW swap the first of pair with Pth

                  KP = ABS( IPIV( K+1 ) )
                  if (KP != K+1) CALL CSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );
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
            if (K < 1) GO TO 60;

            // Test the pivot index.  If greater than zero, a 1 x 1
            // pivot was used, otherwise a 2 x 2 pivot was used.

            if ( IPIV( K ) > 0 ) {

               // 1 x 1 pivot block:

               // Multiply by the diagonal element if forming L * D.

               if (NOUNIT) CALL CSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K != N ) {
                  KP = IPIV( K )

                  // Apply the transformation.

                  cgeru(N-K, NRHS, CONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K - 1

            } else {

               // 2 x 2 pivot block:

               // Multiply by the diagonal block if forming L * D.

               if ( NOUNIT ) {
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D21 = A( K, K-1 )
                  D12 = CONJG( D21 )
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

                  cgeru(N-K, NRHS, CONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );
                  cgeru(N-K, NRHS, CONE, A( K+1, K-1 ), 1, B( K-1, 1 ), LDB, B( K+1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.


                  // Swap the second of pair with IMAXth

                  KP = ABS( IPIV( K ) )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // NOW swap the first of pair with Pth

                  KP = ABS( IPIV( K-1 ) )
                  if (KP != K-1) CALL CSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );

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
   70       IF( K < 1 ) GO TO 90

            // 1 x 1 pivot block.

            if ( IPIV( K ) > 0 ) {
               if ( K > 1 ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation
                     // y = y - B' conjg(x),
                  // where x is a column of A and y is a row of B.

                  clacgv(NRHS, B( K, 1 ), LDB );
                  cgemv('Conjugate', K-1, NRHS, CONE, B, LDB, A( 1, K ), 1, CONE, B( K, 1 ), LDB );
                  clacgv(NRHS, B( K, 1 ), LDB );
               }
               if (NOUNIT) CALL CSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );
               K = K - 1

            // 2 x 2 pivot block.

            } else {
               if ( K > 2 ) {

                  // Swap the second of pair with Pth

                  KP = ABS( IPIV( K ) )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Now swap the first of pair with IMAX(r)th

                  KP = ABS( IPIV( K-1 ) )
                  if (KP != K-1) CALL CSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformations
                     // y = y - B' conjg(x),
                  // where x is a block column of A and y is a block
                  // row of B.

                  clacgv(NRHS, B( K, 1 ), LDB );
                  cgemv('Conjugate', K-2, NRHS, CONE, B, LDB, A( 1, K ), 1, CONE, B( K, 1 ), LDB );
                  clacgv(NRHS, B( K, 1 ), LDB );

                  clacgv(NRHS, B( K-1, 1 ), LDB );
                  cgemv('Conjugate', K-2, NRHS, CONE, B, LDB, A( 1, K-1 ), 1, CONE, B( K-1, 1 ), LDB );
                  clacgv(NRHS, B( K-1, 1 ), LDB );
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D12 = A( K-1, K )
                  D21 = CONJG( D12 )
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
            if (K > N) GO TO 120;

            // 1 x 1 pivot block

            if ( IPIV( K ) > 0 ) {
               if ( K < N ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  clacgv(NRHS, B( K, 1 ), LDB );
                  cgemv('Conjugate', N-K, NRHS, CONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, CONE, B( K, 1 ), LDB );
                  clacgv(NRHS, B( K, 1 ), LDB );
               }
               if (NOUNIT) CALL CSCAL( NRHS, A( K, K ), B( K, 1 ), LDB );
               K = K + 1

            // 2 x 2 pivot block.

            } else {
               if ( K < N-1 ) {

                  // Swap the first of pair with Pth

                  KP = ABS( IPIV( K ) )
                  if (KP != K) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Now swap the second of pair with IMAX(r)th

                  KP = ABS( IPIV( K+1 ) )
                  if (KP != K+1) CALL CSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  clacgv(NRHS, B( K+1, 1 ), LDB );
                  cgemv('Conjugate', N-K-1, NRHS, CONE, B( K+2, 1 ), LDB, A( K+2, K+1 ), 1, CONE, B( K+1, 1 ), LDB );
                  clacgv(NRHS, B( K+1, 1 ), LDB );

                  clacgv(NRHS, B( K, 1 ), LDB );
                  cgemv('Conjugate', N-K-1, NRHS, CONE, B( K+2, 1 ), LDB, A( K+2, K ), 1, CONE, B( K, 1 ), LDB );
                  clacgv(NRHS, B( K, 1 ), LDB );
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D21 = A( K+1, K )
                  D12 = CONJG( D21 )
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

      // End of CLAVHE_ROOK

      }

      SUBROUTINE SLAVSP( UPLO, TRANS, DIAG, N, NRHS, A, IPIV, B, LDB, INFO );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT;
      int                J, K, KC, KCNEXT, KP;
      REAL               D11, D12, D21, D22, T1, T2;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SGER, SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !LSAME( TRANS, 'N' ) && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !LSAME( DIAG, 'U' ) && !LSAME( DIAG, 'N' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('SLAVSP ', -INFO );
         return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOUNIT = LSAME( DIAG, 'N' );
// ------------------------------------------

      // Compute  B := A * B  (No transpose)

// ------------------------------------------
      if ( LSAME( TRANS, 'N' ) ) {

         // Compute  B := U*B
         // where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))

         if ( LSAME( UPLO, 'U' ) ) {

         // Loop forward applying the transformations.

            K = 1;
            KC = 1;
            } // 10
            if (K > N) GO TO 30;

            // 1 x 1 pivot block

            if ( IPIV( K ) > 0 ) {

               // Multiply by the diagonal element if forming U * D.

               if (NOUNIT) CALL SSCAL( NRHS, A( KC+K-1 ), B( K, 1 ), LDB );

               // Multiply by P(K) * inv(U(K))  if K > 1.

               if ( K > 1 ) {

                  // Apply the transformation.

                  sger(K-1, NRHS, ONE, A( KC ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );

                  // Interchange if P(K) != I.

                  KP = IPIV( K );
                  if (KP != K) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               KC = KC + K;
               K = K + 1;
            } else {

               // 2 x 2 pivot block

               KCNEXT = KC + K;

               // Multiply by the diagonal block if forming U * D.

               if ( NOUNIT ) {
                  D11 = A( KCNEXT-1 );
                  D22 = A( KCNEXT+K );
                  D12 = A( KCNEXT+K-1 );
                  D21 = D12;
                  for (J = 1; J <= NRHS; J++) { // 20
                     T1 = B( K, J );
                     T2 = B( K+1, J );
                     B( K, J ) = D11*T1 + D12*T2;
                     B( K+1, J ) = D21*T1 + D22*T2;
                  } // 20
               }

               // Multiply by  P(K) * inv(U(K))  if K > 1.

               if ( K > 1 ) {

                  // Apply the transformations.

                  sger(K-1, NRHS, ONE, A( KC ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB );
                  sger(K-1, NRHS, ONE, A( KCNEXT ), 1, B( K+1, 1 ), LDB, B( 1, 1 ), LDB );

                  // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) );
                  if (KP != K) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               KC = KCNEXT + K + 1;
               K = K + 2;
            }
            GO TO 10;
            } // 30

         // Compute  B := L*B
         // where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .

         } else {

            // Loop backward applying the transformations to B.

            K = N;
            KC = N*( N+1 ) / 2 + 1;
            } // 40
            if (K < 1) GO TO 60;
            KC = KC - ( N-K+1 );

            // Test the pivot index.  If greater than zero, a 1 x 1
            // pivot was used, otherwise a 2 x 2 pivot was used.

            if ( IPIV( K ) > 0 ) {

               // 1 x 1 pivot block:

               // Multiply by the diagonal element if forming L * D.

               if (NOUNIT) CALL SSCAL( NRHS, A( KC ), B( K, 1 ), LDB );

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K != N ) {
                  KP = IPIV( K );

                  // Apply the transformation.

                  sger(N-K, NRHS, ONE, A( KC+1 ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  if (KP != K) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               K = K - 1;

            } else {

               // 2 x 2 pivot block:

               KCNEXT = KC - ( N-K+2 );

               // Multiply by the diagonal block if forming L * D.

               if ( NOUNIT ) {
                  D11 = A( KCNEXT );
                  D22 = A( KC );
                  D21 = A( KCNEXT+1 );
                  D12 = D21;
                  for (J = 1; J <= NRHS; J++) { // 50
                     T1 = B( K-1, J );
                     T2 = B( K, J );
                     B( K-1, J ) = D11*T1 + D12*T2;
                     B( K, J ) = D21*T1 + D22*T2;
                  } // 50
               }

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K != N ) {

                  // Apply the transformation.

                  sger(N-K, NRHS, ONE, A( KC+1 ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB );
                  sger(N-K, NRHS, ONE, A( KCNEXT+2 ), 1, B( K-1, 1 ), LDB, B( K+1, 1 ), LDB );

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  KP = ABS( IPIV( K ) );
                  if (KP != K) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               }
               KC = KCNEXT;
               K = K - 2;
            }
            GO TO 40;
            } // 60
         }
// ----------------------------------------

      // Compute  B := A' * B  (transpose)

// ----------------------------------------
      } else {

         // Form  B := U'*B
         // where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
         // and   U' = inv(U'(1))*P(1)* ... *inv(U'(m))*P(m)

         if ( LSAME( UPLO, 'U' ) ) {

            // Loop backward applying the transformations.

            K = N;
            KC = N*( N+1 ) / 2 + 1;
            } // 70
            if (K < 1) GO TO 90;
            KC = KC - K;

            // 1 x 1 pivot block.

            if ( IPIV( K ) > 0 ) {
               if ( K > 1 ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K );
                  if (KP != K) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  sgemv('Transpose', K-1, NRHS, ONE, B, LDB, A( KC ), 1, ONE, B( K, 1 ), LDB );
               }
               if (NOUNIT) CALL SSCAL( NRHS, A( KC+K-1 ), B( K, 1 ), LDB );
               K = K - 1;

            // 2 x 2 pivot block.

            } else {
               KCNEXT = KC - ( K-1 );
               if ( K > 2 ) {

                  // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) );
                  if (KP != K-1) CALL SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformations

                  sgemv('Transpose', K-2, NRHS, ONE, B, LDB, A( KC ), 1, ONE, B( K, 1 ), LDB );
                  sgemv('Transpose', K-2, NRHS, ONE, B, LDB, A( KCNEXT ), 1, ONE, B( K-1, 1 ), LDB );
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( KC-1 );
                  D22 = A( KC+K-1 );
                  D12 = A( KC+K-2 );
                  D21 = D12;
                  for (J = 1; J <= NRHS; J++) { // 80
                     T1 = B( K-1, J );
                     T2 = B( K, J );
                     B( K-1, J ) = D11*T1 + D12*T2;
                     B( K, J ) = D21*T1 + D22*T2;
                  } // 80
               }
               KC = KCNEXT;
               K = K - 2;
            }
            GO TO 70;
            } // 90

         // Form  B := L'*B
         // where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
         // and   L' = inv(L(m))*P(m)* ... *inv(L(1))*P(1)

         } else {

            // Loop forward applying the L-transformations.

            K = 1;
            KC = 1;
            } // 100
            if (K > N) GO TO 120;

            // 1 x 1 pivot block

            if ( IPIV( K ) > 0 ) {
               if ( K < N ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K );
                  if (KP != K) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  sgemv('Transpose', N-K, NRHS, ONE, B( K+1, 1 ), LDB, A( KC+1 ), 1, ONE, B( K, 1 ), LDB );
               }
               if (NOUNIT) CALL SSCAL( NRHS, A( KC ), B( K, 1 ), LDB );
               KC = KC + N - K + 1;
               K = K + 1;

            // 2 x 2 pivot block.

            } else {
               KCNEXT = KC + N - K + 1;
               if ( K < N-1 ) {

               // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) );
                  if (KP != K+1) CALL SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB );

                  // Apply the transformation

                  sgemv('Transpose', N-K-1, NRHS, ONE, B( K+2, 1 ), LDB, A( KCNEXT+1 ), 1, ONE, B( K+1, 1 ), LDB );
                  sgemv('Transpose', N-K-1, NRHS, ONE, B( K+2, 1 ), LDB, A( KC+2 ), 1, ONE, B( K, 1 ), LDB );
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( KC );
                  D22 = A( KCNEXT );
                  D21 = A( KC+1 );
                  D12 = D21;
                  for (J = 1; J <= NRHS; J++) { // 110
                     T1 = B( K, J );
                     T2 = B( K+1, J );
                     B( K, J ) = D11*T1 + D12*T2;
                     B( K+1, J ) = D21*T1 + D22*T2;
                  } // 110
               }
               KC = KCNEXT + ( N-K );
               K = K + 2;
            }
            GO TO 100;
            } // 120
         }

      }
      return;

      // End of SLAVSP

      }

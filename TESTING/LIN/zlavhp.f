      SUBROUTINE ZLAVHP( UPLO, TRANS, DIAG, N, NRHS, A, IPIV, B, LDB, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT;
      int                J, K, KC, KCNEXT, KP;
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
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.LSAME( DIAG, 'U' ) .AND. .NOT.LSAME( DIAG, 'N' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZLAVHP ', -INFO )
         RETURN
      }

      // Quick return if possible.

      IF( N.EQ.0 ) RETURN

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
            KC = 1
   10       CONTINUE
            IF( K.GT.N ) GO TO 30

            // 1 x 1 pivot block

            if ( IPIV( K ).GT.0 ) {

               // Multiply by the diagonal element if forming U * D.

               IF( NOUNIT ) CALL ZSCAL( NRHS, A( KC+K-1 ), B( K, 1 ), LDB )

               // Multiply by P(K) * inv(U(K))  if K > 1.

               if ( K.GT.1 ) {

                  // Apply the transformation.

                  CALL ZGERU( K-1, NRHS, ONE, A( KC ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB )

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               }
               KC = KC + K
               K = K + 1
            } else {

               // 2 x 2 pivot block

               KCNEXT = KC + K

               // Multiply by the diagonal block if forming U * D.

               if ( NOUNIT ) {
                  D11 = A( KCNEXT-1 )
                  D22 = A( KCNEXT+K )
                  D12 = A( KCNEXT+K-1 )
                  D21 = DCONJG( D12 )
                  DO 20 J = 1, NRHS
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
   20             CONTINUE
               }

               // Multiply by  P(K) * inv(U(K))  if K > 1.

               if ( K.GT.1 ) {

                  // Apply the transformations.

                  CALL ZGERU( K-1, NRHS, ONE, A( KC ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB )                   CALL ZGERU( K-1, NRHS, ONE, A( KCNEXT ), 1, B( K+1, 1 ), LDB, B( 1, 1 ), LDB )

                  // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               }
               KC = KCNEXT + K + 1
               K = K + 2
            }
            GO TO 10
   30       CONTINUE

         // Compute  B := L*B
         // where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .

         } else {

            // Loop backward applying the transformations to B.

            K = N
            KC = N*( N+1 ) / 2 + 1
   40       CONTINUE
            IF( K.LT.1 ) GO TO 60
            KC = KC - ( N-K+1 )

            // Test the pivot index.  If greater than zero, a 1 x 1
            // pivot was used, otherwise a 2 x 2 pivot was used.

            if ( IPIV( K ).GT.0 ) {

               // 1 x 1 pivot block:

               // Multiply by the diagonal element if forming L * D.

               IF( NOUNIT ) CALL ZSCAL( NRHS, A( KC ), B( K, 1 ), LDB )

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K.NE.N ) {
                  KP = IPIV( K )

                  // Apply the transformation.

                  CALL ZGERU( N-K, NRHS, ONE, A( KC+1 ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB )

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               }
               K = K - 1

            } else {

               // 2 x 2 pivot block:

               KCNEXT = KC - ( N-K+2 )

               // Multiply by the diagonal block if forming L * D.

               if ( NOUNIT ) {
                  D11 = A( KCNEXT )
                  D22 = A( KC )
                  D21 = A( KCNEXT+1 )
                  D12 = DCONJG( D21 )
                  DO 50 J = 1, NRHS
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
   50             CONTINUE
               }

               // Multiply by  P(K) * inv(L(K))  if K < N.

               if ( K.NE.N ) {

                  // Apply the transformation.

                  CALL ZGERU( N-K, NRHS, ONE, A( KC+1 ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB )                   CALL ZGERU( N-K, NRHS, ONE, A( KCNEXT+2 ), 1, B( K-1, 1 ), LDB, B( K+1, 1 ), LDB )

                  // Interchange if a permutation was applied at the
                  // K-th step of the factorization.

                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               }
               KC = KCNEXT
               K = K - 2
            }
            GO TO 40
   60       CONTINUE
         }
*-------------------------------------------------

      // Compute  B := A^H * B  (conjugate transpose)

*-------------------------------------------------
      } else {

         // Form  B := U^H*B
         // where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
         // and   U^H = inv(U^H(1))*P(1)* ... *inv(U^H(m))*P(m)

         if ( LSAME( UPLO, 'U' ) ) {

            // Loop backward applying the transformations.

            K = N
            KC = N*( N+1 ) / 2 + 1
   70       CONTINUE
            IF( K.LT.1 ) GO TO 90
            KC = KC - K

            // 1 x 1 pivot block.

            if ( IPIV( K ).GT.0 ) {
               if ( K.GT.1 ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )

                  // Apply the transformation:
                     // y := y - B' * conjg(x)
                  // where x is a column of A and y is a row of B.

                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-1, NRHS, ONE, B, LDB, A( KC ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               }
               IF( NOUNIT ) CALL ZSCAL( NRHS, A( KC+K-1 ), B( K, 1 ), LDB )
               K = K - 1

            // 2 x 2 pivot block.

            } else {
               KCNEXT = KC - ( K-1 )
               if ( K.GT.2 ) {

                  // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K-1 ) CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )

                  // Apply the transformations.

                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-2, NRHS, ONE, B, LDB, A( KC ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )

                  CALL ZLACGV( NRHS, B( K-1, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-2, NRHS, ONE, B, LDB, A( KCNEXT ), 1, ONE, B( K-1, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K-1, 1 ), LDB )
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( KC-1 )
                  D22 = A( KC+K-1 )
                  D12 = A( KC+K-2 )
                  D21 = DCONJG( D12 )
                  DO 80 J = 1, NRHS
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
   80             CONTINUE
               }
               KC = KCNEXT
               K = K - 2
            }
            GO TO 70
   90       CONTINUE

         // Form  B := L^H*B
         // where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
         // and   L^H = inv(L(m))*P(m)* ... *inv(L(1))*P(1)

         } else {

            // Loop forward applying the L-transformations.

            K = 1
            KC = 1
  100       CONTINUE
            IF( K.GT.N ) GO TO 120

            // 1 x 1 pivot block

            if ( IPIV( K ).GT.0 ) {
               if ( K.LT.N ) {

                  // Interchange if P(K) != I.

                  KP = IPIV( K )
                  IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )

                  // Apply the transformation

                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K, NRHS, ONE, B( K+1, 1 ), LDB, A( KC+1 ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               }
               IF( NOUNIT ) CALL ZSCAL( NRHS, A( KC ), B( K, 1 ), LDB )
               KC = KC + N - K + 1
               K = K + 1

            // 2 x 2 pivot block.

            } else {
               KCNEXT = KC + N - K + 1
               if ( K.LT.N-1 ) {

               // Interchange if P(K) != I.

                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K+1 ) CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )

                  // Apply the transformation

                  CALL ZLACGV( NRHS, B( K+1, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K-1, NRHS, ONE, B( K+2, 1 ), LDB, A( KCNEXT+1 ), 1, ONE, B( K+1, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K+1, 1 ), LDB )

                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K-1, NRHS, ONE, B( K+2, 1 ), LDB, A( KC+2 ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               }

               // Multiply by the diagonal block if non-unit.

               if ( NOUNIT ) {
                  D11 = A( KC )
                  D22 = A( KCNEXT )
                  D21 = A( KC+1 )
                  D12 = DCONJG( D21 )
                  DO 110 J = 1, NRHS
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
  110             CONTINUE
               }
               KC = KCNEXT + ( N-K )
               K = K + 2
            }
            GO TO 100
  120       CONTINUE
         }

      }
      RETURN

      // End of ZLAVHP

      }

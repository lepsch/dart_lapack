      SUBROUTINE ZLAVHE_ROOK( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV,
     $                        B, LDB, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
      INTEGER            J, K, KP
      COMPLEX*16         D11, D12, D21, D22, T1, T2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERU, ZLACGV, ZSCAL, ZSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) )
     $          THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( DIAG, 'U' ) .AND. .NOT.LSAME( DIAG, 'N' ) )
     $          THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAVHE_ROOK ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*------------------------------------------
*
*     Compute  B := A * B  (No transpose)
*
*------------------------------------------
      IF( LSAME( TRANS, 'N' ) ) THEN
*
*        Compute  B := U*B
*        where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Loop forward applying the transformations.
*
            K = 1
   10       CONTINUE
            IF( K.GT.N )
     $         GO TO 30
            IF( IPIV( K ).GT.0 ) THEN
*
*              1 x 1 pivot block
*
*              Multiply by the diagonal element if forming U * D.
*
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
*
*              Multiply by  P(K) * inv(U(K))  if K > 1.
*
               IF( K.GT.1 ) THEN
*
*                 Apply the transformation.
*
                  CALL ZGERU( K-1, NRHS, CONE, A( 1, K ), 1, B( K, 1 ),
     $                        LDB, B( 1, 1 ), LDB )
*
*                 Interchange if P(K) != I.
*
                  KP = IPIV( K )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               END IF
               K = K + 1
            ELSE
*
*              2 x 2 pivot block
*
*              Multiply by the diagonal block if forming U * D.
*
               IF( NOUNIT ) THEN
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D12 = A( K, K+1 )
                  D21 = DCONJG( D12 )
                  DO 20 J = 1, NRHS
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
   20             CONTINUE
               END IF
*
*              Multiply by  P(K) * inv(U(K))  if K > 1.
*
               IF( K.GT.1 ) THEN
*
*                 Apply the transformations.
*
                  CALL ZGERU( K-1, NRHS, CONE, A( 1, K ), 1, B( K, 1 ),
     $                        LDB, B( 1, 1 ), LDB )
                  CALL ZGERU( K-1, NRHS, CONE, A( 1, K+1 ), 1,
     $                        B( K+1, 1 ), LDB, B( 1, 1 ), LDB )
*
*                 Interchange if a permutation was applied at the
*                 K-th step of the factorization.
*
*                 Swap the first of pair with IMAXth
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 NOW swap the first of pair with Pth
*
                  KP = ABS( IPIV( K+1 ) )
                  IF( KP.NE.K+1 )
     $               CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ),
     $                           LDB )
               END IF
               K = K + 2
            END IF
            GO TO 10
   30       CONTINUE
*
*        Compute  B := L*B
*        where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
*
         ELSE
*
*           Loop backward applying the transformations to B.
*
            K = N
   40       CONTINUE
            IF( K.LT.1 )
     $         GO TO 60
*
*           Test the pivot index.  If greater than zero, a 1 x 1
*           pivot was used, otherwise a 2 x 2 pivot was used.
*
            IF( IPIV( K ).GT.0 ) THEN
*
*              1 x 1 pivot block:
*
*              Multiply by the diagonal element if forming L * D.
*
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
*
*              Multiply by  P(K) * inv(L(K))  if K < N.
*
               IF( K.NE.N ) THEN
                  KP = IPIV( K )
*
*                 Apply the transformation.
*
                  CALL ZGERU( N-K, NRHS, CONE, A( K+1, K ), 1,
     $                        B( K, 1 ), LDB, B( K+1, 1 ), LDB )
*
*                 Interchange if a permutation was applied at the
*                 K-th step of the factorization.
*
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               END IF
               K = K - 1
*
            ELSE
*
*              2 x 2 pivot block:
*
*              Multiply by the diagonal block if forming L * D.
*
               IF( NOUNIT ) THEN
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D21 = A( K, K-1 )
                  D12 = DCONJG( D21 )
                  DO 50 J = 1, NRHS
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
   50             CONTINUE
               END IF
*
*              Multiply by  P(K) * inv(L(K))  if K < N.
*
               IF( K.NE.N ) THEN
*
*                 Apply the transformation.
*
                  CALL ZGERU( N-K, NRHS, CONE, A( K+1, K ), 1,
     $                        B( K, 1 ), LDB, B( K+1, 1 ), LDB )
                  CALL ZGERU( N-K, NRHS, CONE, A( K+1, K-1 ), 1,
     $                        B( K-1, 1 ), LDB, B( K+1, 1 ), LDB )
*
*                 Interchange if a permutation was applied at the
*                 K-th step of the factorization.
*
*
*                 Swap the second of pair with IMAXth
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 NOW swap the first of pair with Pth
*
                  KP = ABS( IPIV( K-1 ) )
                  IF( KP.NE.K-1 )
     $               CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ),
     $                           LDB )
*
               END IF
               K = K - 2
            END IF
            GO TO 40
   60       CONTINUE
         END IF
*--------------------------------------------------
*
*     Compute  B := A^H * B  (conjugate transpose)
*
*--------------------------------------------------
      ELSE
*
*        Form  B := U^H*B
*        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
*        and   U^H = inv(U^H(1))*P(1)* ... *inv(U^H(m))*P(m)
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Loop backward applying the transformations.
*
            K = N
   70       IF( K.LT.1 )
     $         GO TO 90
*
*           1 x 1 pivot block.
*
            IF( IPIV( K ).GT.0 ) THEN
               IF( K.GT.1 ) THEN
*
*                 Interchange if P(K) != I.
*
                  KP = IPIV( K )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 Apply the transformation
*                    y = y - B' DCONJG(x),
*                 where x is a column of A and y is a row of B.
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-1, NRHS, CONE, B, LDB,
     $                        A( 1, K ), 1, CONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               END IF
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
               K = K - 1
*
*           2 x 2 pivot block.
*
            ELSE
               IF( K.GT.2 ) THEN
*
*                 Swap the second of pair with Pth
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 Now swap the first of pair with IMAX(r)th
*
                  KP = ABS( IPIV( K-1 ) )
                  IF( KP.NE.K-1 )
     $               CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ),
     $                           LDB )
*
*                 Apply the transformations
*                    y = y - B' DCONJG(x),
*                 where x is a block column of A and y is a block
*                 row of B.
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-2, NRHS, CONE, B, LDB,
     $                        A( 1, K ), 1, CONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
*
                  CALL ZLACGV( NRHS, B( K-1, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-2, NRHS, CONE, B, LDB,
     $                        A( 1, K-1 ), 1, CONE, B( K-1, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K-1, 1 ), LDB )
               END IF
*
*              Multiply by the diagonal block if non-unit.
*
               IF( NOUNIT ) THEN
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D12 = A( K-1, K )
                  D21 = DCONJG( D12 )
                  DO 80 J = 1, NRHS
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
   80             CONTINUE
               END IF
               K = K - 2
            END IF
            GO TO 70
   90       CONTINUE
*
*        Form  B := L^H*B
*        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
*        and   L^H = inv(L^H(m))*P(m)* ... *inv(L^H(1))*P(1)
*
         ELSE
*
*           Loop forward applying the L-transformations.
*
            K = 1
  100       CONTINUE
            IF( K.GT.N )
     $         GO TO 120
*
*           1 x 1 pivot block
*
            IF( IPIV( K ).GT.0 ) THEN
               IF( K.LT.N ) THEN
*
*                 Interchange if P(K) != I.
*
                  KP = IPIV( K )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 Apply the transformation
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K, NRHS, CONE, B( K+1, 1 ),
     $                       LDB, A( K+1, K ), 1, CONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               END IF
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
               K = K + 1
*
*           2 x 2 pivot block.
*
            ELSE
               IF( K.LT.N-1 ) THEN
*
*                 Swap the first of pair with Pth
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 Now swap the second of pair with IMAX(r)th
*
                  KP = ABS( IPIV( K+1 ) )
                  IF( KP.NE.K+1 )
     $               CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ),
     $                           LDB )
*
*                 Apply the transformation
*
                  CALL ZLACGV( NRHS, B( K+1, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K-1, NRHS, CONE,
     $                        B( K+2, 1 ), LDB, A( K+2, K+1 ), 1, CONE,
     $                        B( K+1, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K+1, 1 ), LDB )
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K-1, NRHS, CONE,
     $                        B( K+2, 1 ), LDB, A( K+2, K ), 1, CONE,
     $                        B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               END IF
*
*              Multiply by the diagonal block if non-unit.
*
               IF( NOUNIT ) THEN
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D21 = A( K+1, K )
                  D12 = DCONJG( D21 )
                  DO 110 J = 1, NRHS
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
  110             CONTINUE
               END IF
               K = K + 2
            END IF
            GO TO 100
  120       CONTINUE
         END IF
*
      END IF
      RETURN
*
*     End of ZLAVHE_ROOK
*
      END
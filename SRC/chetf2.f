      SUBROUTINE CHETF2( UPLO, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX            A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0E+0, SEVTEN = 17.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      int                I, IMAX, J, JMAX, K, KK, KP, KSTEP
      REAL               ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT
      COMPLEX            D12, D21, T, WK, WKM1, WKP1, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      int                ICAMAX
      REAL               SLAPY2
      EXTERNAL           LSAME, ICAMAX, SLAPY2, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHER, CSSCAL, CSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, MAX, REAL, SQRT
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHETF2', -INFO )
         RETURN
      END IF
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U**H using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 ) GO TO 90
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( REAL( A( K, K ) ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.GT.1 ) THEN
            IMAX = ICAMAX( K-1, A( 1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN
*
*           Column K is or underflow, or contains a NaN:
*           set INFO and continue
*
            IF( INFO.EQ.0 ) INFO = K
            KP = K
            A( K, K ) = REAL( A( K, K ) )
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = IMAX + ICAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = CABS1( A( IMAX, JMAX ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = ICAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( ABS( REAL( A( IMAX, IMAX ) ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
               ELSE
*
*                 interchange rows and columns K-1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the leading
*              submatrix A(1:k,1:k)
*
               CALL CSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               DO 20 J = KP + 1, KK - 1
                  T = CONJG( A( J, KK ) )
                  A( J, KK ) = CONJG( A( KP, J ) )
                  A( KP, J ) = T
   20          CONTINUE
               A( KP, KK ) = CONJG( A( KP, KK ) )
               R1 = REAL( A( KK, KK ) )
               A( KK, KK ) = REAL( A( KP, KP ) )
               A( KP, KP ) = R1
               IF( KSTEP.EQ.2 ) THEN
                  A( K, K ) = REAL( A( K, K ) )
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            ELSE
               A( K, K ) = REAL( A( K, K ) )
               IF( KSTEP.EQ.2 ) A( K-1, K-1 ) = REAL( A( K-1, K-1 ) )
            END IF
*
*           Update the leading submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
*              Perform a rank-1 update of A(1:k-1,1:k-1) as
*
*              A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H
*
               R1 = ONE / REAL( A( K, K ) )
               CALL CHER( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
*
*              Store U(k) in column k
*
               CALL CSSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k-1 now hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
*              Perform a rank-2 update of A(1:k-2,1:k-2) as
*
*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H
*
               IF( K.GT.2 ) THEN
*
                  D = SLAPY2( REAL( A( K-1, K ) ), AIMAG( A( K-1, K ) ) )
                  D22 = REAL( A( K-1, K-1 ) ) / D
                  D11 = REAL( A( K, K ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D12 = A( K-1, K ) / D
                  D = TT / D
*
                  DO 40 J = K - 2, 1, -1
                     WKM1 = D*( D11*A( J, K-1 )-CONJG( D12 )*A( J, K ) )
                     WK = D*( D22*A( J, K )-D12*A( J, K-1 ) )
                     DO 30 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*CONJG( WK ) - A( I, K-1 )*CONJG( WKM1 )
   30                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
                     A( J, J ) = CMPLX( REAL( A( J, J ) ), 0.0E+0 )
   40             CONTINUE
*
               END IF
*
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
      ELSE
*
*        Factorize A as L*D*L**H using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2
*
         K = 1
   50    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N ) GO TO 90
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = ABS( REAL( A( K, K ) ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.LT.N ) THEN
            IMAX = K + ICAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN
*
*           Column K is zero or underflow, contains a NaN:
*           set INFO and continue
*
            IF( INFO.EQ.0 ) INFO = K
            KP = K
            A( K, K ) = REAL( A( K, K ) )
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = K - 1 + ICAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = CABS1( A( IMAX, JMAX ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + ICAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( ABS( REAL( A( IMAX, IMAX ) ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
               ELSE
*
*                 interchange rows and columns K+1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the trailing
*              submatrix A(k:n,k:n)
*
               IF( KP.LT.N ) CALL CSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               DO 60 J = KK + 1, KP - 1
                  T = CONJG( A( J, KK ) )
                  A( J, KK ) = CONJG( A( KP, J ) )
                  A( KP, J ) = T
   60          CONTINUE
               A( KP, KK ) = CONJG( A( KP, KK ) )
               R1 = REAL( A( KK, KK ) )
               A( KK, KK ) = REAL( A( KP, KP ) )
               A( KP, KP ) = R1
               IF( KSTEP.EQ.2 ) THEN
                  A( K, K ) = REAL( A( K, K ) )
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            ELSE
               A( K, K ) = REAL( A( K, K ) )
               IF( KSTEP.EQ.2 ) A( K+1, K+1 ) = REAL( A( K+1, K+1 ) )
            END IF
*
*           Update the trailing submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
               IF( K.LT.N ) THEN
*
*                 Perform a rank-1 update of A(k+1:n,k+1:n) as
*
*                 A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H
*
                  R1 = ONE / REAL( A( K, K ) )
                  CALL CHER( UPLO, N-K, -R1, A( K+1, K ), 1, A( K+1, K+1 ), LDA )
*
*                 Store L(k) in column K
*
                  CALL CSSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-by-2 pivot block D(k)
*
               IF( K.LT.N-1 ) THEN
*
*                 Perform a rank-2 update of A(k+2:n,k+2:n) as
*
*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H
*
*                 where L(k) and L(k+1) are the k-th and (k+1)-th
*                 columns of L
*
                  D = SLAPY2( REAL( A( K+1, K ) ), AIMAG( A( K+1, K ) ) )
                  D11 = REAL( A( K+1, K+1 ) ) / D
                  D22 = REAL( A( K, K ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D21 = A( K+1, K ) / D
                  D =  TT / D
*
                  DO 80 J = K + 2, N
                     WK = D*( D11*A( J, K )-D21*A( J, K+1 ) )
                     WKP1 = D*( D22*A( J, K+1 )-CONJG( D21 )*A( J, K ) )
                     DO 70 I = J, N
                        A( I, J ) = A( I, J ) - A( I, K )*CONJG( WK ) - A( I, K+1 )*CONJG( WKP1 )
   70                CONTINUE
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
                     A( J, J ) = CMPLX( REAL( A( J, J ) ), 0.0E+0 )
   80             CONTINUE
               END IF
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 50
*
      END IF
*
   90 CONTINUE
      RETURN
*
*     End of CHETF2
*
      END

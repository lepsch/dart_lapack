      SUBROUTINE CHPTRF( UPLO, N, AP, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      REAL               EIGHT, SEVTEN
      const              EIGHT = 8.0E+0, SEVTEN = 17.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IMAX, J, JMAX, K, KC, KK, KNC, KP, KPC, KSTEP, KX, NPP;
      REAL               ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT;
      COMPLEX            D12, D21, T, WK, WKM1, WKP1, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAPY2
      // EXTERNAL LSAME, ICAMAX, SLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPR, CSSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHPTRF', -INFO )
         RETURN
      END IF

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      IF( UPPER ) THEN

         // Factorize A as U*D*U**H using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N
         KC = ( N-1 )*N / 2 + 1
   10    CONTINUE
         KNC = KC

         // If K < 1, exit from loop

         IF( K.LT.1 ) GO TO 110
         KSTEP = 1

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( REAL( AP( KC+K-1 ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value

         IF( K.GT.1 ) THEN
            IMAX = ICAMAX( K-1, AP( KC ), 1 )
            COLMAX = CABS1( AP( KC+IMAX-1 ) )
         ELSE
            COLMAX = ZERO
         END IF

         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN

            // Column K is zero: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            AP( KC+K-1 ) = REAL( AP( KC+K-1 ) )
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN

               // no interchange, use 1-by-1 pivot block

               KP = K
            ELSE

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               ROWMAX = ZERO
               JMAX = IMAX
               KX = IMAX*( IMAX+1 ) / 2 + IMAX
               DO 20 J = IMAX + 1, K
                  IF( CABS1( AP( KX ) ).GT.ROWMAX ) THEN
                     ROWMAX = CABS1( AP( KX ) )
                     JMAX = J
                  END IF
                  KX = KX + J
   20          CONTINUE
               KPC = ( IMAX-1 )*IMAX / 2 + 1
               IF( IMAX.GT.1 ) THEN
                  JMAX = ICAMAX( IMAX-1, AP( KPC ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( AP( KPC+JMAX-1 ) ) )
               END IF

               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN

                  // no interchange, use 1-by-1 pivot block

                  KP = K
               ELSE IF( ABS( REAL( AP( KPC+IMAX-1 ) ) ).GE.ALPHA* ROWMAX ) THEN

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX
               ELSE

                  // interchange rows and columns K-1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF

            KK = K - KSTEP + 1
            IF( KSTEP.EQ.2 ) KNC = KNC - K + 1
            IF( KP.NE.KK ) THEN

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               CALL CSWAP( KP-1, AP( KNC ), 1, AP( KPC ), 1 )
               KX = KPC + KP - 1
               DO 30 J = KP + 1, KK - 1
                  KX = KX + J - 1
                  T = CONJG( AP( KNC+J-1 ) )
                  AP( KNC+J-1 ) = CONJG( AP( KX ) )
                  AP( KX ) = T
   30          CONTINUE
               AP( KX+KK-1 ) = CONJG( AP( KX+KK-1 ) )
               R1 = REAL( AP( KNC+KK-1 ) )
               AP( KNC+KK-1 ) = REAL( AP( KPC+KP-1 ) )
               AP( KPC+KP-1 ) = R1
               IF( KSTEP.EQ.2 ) THEN
                  AP( KC+K-1 ) = REAL( AP( KC+K-1 ) )
                  T = AP( KC+K-2 )
                  AP( KC+K-2 ) = AP( KC+KP-1 )
                  AP( KC+KP-1 ) = T
               END IF
            ELSE
               AP( KC+K-1 ) = REAL( AP( KC+K-1 ) )
               IF( KSTEP.EQ.2 ) AP( KC-1 ) = REAL( AP( KC-1 ) )
            END IF

            // Update the leading submatrix

            IF( KSTEP.EQ.1 ) THEN

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Perform a rank-1 update of A(1:k-1,1:k-1) as

               // A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H

               R1 = ONE / REAL( AP( KC+K-1 ) )
               CALL CHPR( UPLO, K-1, -R1, AP( KC ), 1, AP )

               // Store U(k) in column k

               CALL CSSCAL( K-1, R1, AP( KC ), 1 )
            ELSE

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
                  // = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H

               IF( K.GT.2 ) THEN

                  D = SLAPY2( REAL( AP( K-1+( K-1 )*K / 2 ) ), AIMAG( AP( K-1+( K-1 )*K / 2 ) ) )
                  D22 = REAL( AP( K-1+( K-2 )*( K-1 ) / 2 ) ) / D
                  D11 = REAL( AP( K+( K-1 )*K / 2 ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D12 = AP( K-1+( K-1 )*K / 2 ) / D
                  D = TT / D

                  DO 50 J = K - 2, 1, -1
                     WKM1 = D*( D11*AP( J+( K-2 )*( K-1 ) / 2 )- CONJG( D12 )*AP( J+( K-1 )*K / 2 ) )                      WK = D*( D22*AP( J+( K-1 )*K / 2 )-D12* AP( J+( K-2 )*( K-1 ) / 2 ) )
                     DO 40 I = J, 1, -1
                        AP( I+( J-1 )*J / 2 ) = AP( I+( J-1 )*J / 2 ) - AP( I+( K-1 )*K / 2 )*CONJG( WK ) - AP( I+( K-2 )*( K-1 ) / 2 )*CONJG( WKM1 )
   40                CONTINUE
                     AP( J+( K-1 )*K / 2 ) = WK
                     AP( J+( K-2 )*( K-1 ) / 2 ) = WKM1
                     AP( J+( J-1 )*J / 2 ) = CMPLX( REAL( AP( J+( J-1 )* J / 2 ) ), 0.0E+0 )
   50             CONTINUE

               END IF

            END IF
         END IF

         // Store details of the interchanges in IPIV

         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF

         // Decrease K and return to the start of the main loop

         K = K - KSTEP
         KC = KNC - K
         GO TO 10

      ELSE

         // Factorize A as L*D*L**H using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1
         KC = 1
         NPP = N*( N+1 ) / 2
   60    CONTINUE
         KNC = KC

         // If K > N, exit from loop

         IF( K.GT.N ) GO TO 110
         KSTEP = 1

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( REAL( AP( KC ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value

         IF( K.LT.N ) THEN
            IMAX = K + ICAMAX( N-K, AP( KC+1 ), 1 )
            COLMAX = CABS1( AP( KC+IMAX-K ) )
         ELSE
            COLMAX = ZERO
         END IF

         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN

            // Column K is zero: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            AP( KC ) = REAL( AP( KC ) )
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN

               // no interchange, use 1-by-1 pivot block

               KP = K
            ELSE

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               ROWMAX = ZERO
               KX = KC + IMAX - K
               DO 70 J = K, IMAX - 1
                  IF( CABS1( AP( KX ) ).GT.ROWMAX ) THEN
                     ROWMAX = CABS1( AP( KX ) )
                     JMAX = J
                  END IF
                  KX = KX + N - J
   70          CONTINUE
               KPC = NPP - ( N-IMAX+1 )*( N-IMAX+2 ) / 2 + 1
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + ICAMAX( N-IMAX, AP( KPC+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( AP( KPC+JMAX-IMAX ) ) )
               END IF

               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN

                  // no interchange, use 1-by-1 pivot block

                  KP = K
               ELSE IF( ABS( REAL( AP( KPC ) ) ).GE.ALPHA*ROWMAX ) THEN

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX
               ELSE

                  // interchange rows and columns K+1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF

            KK = K + KSTEP - 1
            IF( KSTEP.EQ.2 ) KNC = KNC + N - K + 1
            IF( KP.NE.KK ) THEN

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               IF( KP.LT.N ) CALL CSWAP( N-KP, AP( KNC+KP-KK+1 ), 1, AP( KPC+1 ), 1 )
               KX = KNC + KP - KK
               DO 80 J = KK + 1, KP - 1
                  KX = KX + N - J + 1
                  T = CONJG( AP( KNC+J-KK ) )
                  AP( KNC+J-KK ) = CONJG( AP( KX ) )
                  AP( KX ) = T
   80          CONTINUE
               AP( KNC+KP-KK ) = CONJG( AP( KNC+KP-KK ) )
               R1 = REAL( AP( KNC ) )
               AP( KNC ) = REAL( AP( KPC ) )
               AP( KPC ) = R1
               IF( KSTEP.EQ.2 ) THEN
                  AP( KC ) = REAL( AP( KC ) )
                  T = AP( KC+1 )
                  AP( KC+1 ) = AP( KC+KP-K )
                  AP( KC+KP-K ) = T
               END IF
            ELSE
               AP( KC ) = REAL( AP( KC ) )
               IF( KSTEP.EQ.2 ) AP( KNC ) = REAL( AP( KNC ) )
            END IF

            // Update the trailing submatrix

            IF( KSTEP.EQ.1 ) THEN

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               IF( K.LT.N ) THEN

                  // Perform a rank-1 update of A(k+1:n,k+1:n) as

                  // A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H

                  R1 = ONE / REAL( AP( KC ) )
                  CALL CHPR( UPLO, N-K, -R1, AP( KC+1 ), 1, AP( KC+N-K+1 ) )

                  // Store L(k) in column K

                  CALL CSSCAL( N-K, R1, AP( KC+1 ), 1 )
               END IF
            ELSE

               // 2-by-2 pivot block D(k): columns K and K+1 now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               IF( K.LT.N-1 ) THEN

                  // Perform a rank-2 update of A(k+2:n,k+2:n) as

                  // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
                     // = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H

                  // where L(k) and L(k+1) are the k-th and (k+1)-th
                  // columns of L

                  D = SLAPY2( REAL( AP( K+1+( K-1 )*( 2*N-K ) / 2 ) ), AIMAG( AP( K+1+( K-1 )*( 2*N-K ) / 2 ) ) )
                  D11 = REAL( AP( K+1+K*( 2*N-K-1 ) / 2 ) ) / D
                  D22 = REAL( AP( K+( K-1 )*( 2*N-K ) / 2 ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D21 = AP( K+1+( K-1 )*( 2*N-K ) / 2 ) / D
                  D = TT / D

                  DO 100 J = K + 2, N
                     WK = D*( D11*AP( J+( K-1 )*( 2*N-K ) / 2 )-D21* AP( J+K*( 2*N-K-1 ) / 2 ) )                      WKP1 = D*( D22*AP( J+K*( 2*N-K-1 ) / 2 )- CONJG( D21 )*AP( J+( K-1 )*( 2*N-K ) / 2 ) )
                     DO 90 I = J, N
                        AP( I+( J-1 )*( 2*N-J ) / 2 ) = AP( I+( J-1 )* ( 2*N-J ) / 2 ) - AP( I+( K-1 )*( 2*N-K ) / 2 )*CONJG( WK ) - AP( I+K*( 2*N-K-1 ) / 2 )* CONJG( WKP1 )
   90                CONTINUE
                     AP( J+( K-1 )*( 2*N-K ) / 2 ) = WK
                     AP( J+K*( 2*N-K-1 ) / 2 ) = WKP1
                     AP( J+( J-1 )*( 2*N-J ) / 2 ) = CMPLX( REAL( AP( J+( J-1 )*( 2*N-J ) / 2 ) ), 0.0E+0 )
  100             CONTINUE
               END IF
            END IF
         END IF

         // Store details of the interchanges in IPIV

         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF

         // Increase K and return to the start of the main loop

         K = K + KSTEP
         KC = KNC + N - K + 2
         GO TO 60

      END IF

  110 CONTINUE
      RETURN

      // End of CHPTRF

      }

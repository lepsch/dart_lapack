      SUBROUTINE ZHPTRF( UPLO, N, AP, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0D+0, SEVTEN = 17.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IMAX, J, JMAX, K, KC, KK, KNC, KP, KPC, KSTEP, KX, NPP;
      double             ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT;
      COMPLEX*16         D12, D21, T, WK, WKM1, WKP1, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAPY2;
      // EXTERNAL LSAME, IZAMAX, DLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZHPR, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('ZHPTRF', -INFO );
         RETURN
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      if ( UPPER ) {

         // Factorize A as U*D*U**H using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N
         KC = ( N-1 )*N / 2 + 1
         } // 10
         KNC = KC

         // If K < 1, exit from loop

         IF( K.LT.1 ) GO TO 110
         KSTEP = 1

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( DBLE( AP( KC+K-1 ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value

         if ( K.GT.1 ) {
            IMAX = IZAMAX( K-1, AP( KC ), 1 )
            COLMAX = CABS1( AP( KC+IMAX-1 ) )
         } else {
            COLMAX = ZERO
         }

         if ( MAX( ABSAKK, COLMAX ).EQ.ZERO ) {

            // Column K is zero: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            AP( KC+K-1 ) = DBLE( AP( KC+K-1 ) )
         } else {
            if ( ABSAKK.GE.ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               ROWMAX = ZERO
               JMAX = IMAX
               KX = IMAX*( IMAX+1 ) / 2 + IMAX
               DO 20 J = IMAX + 1, K
                  if ( CABS1( AP( KX ) ).GT.ROWMAX ) {
                     ROWMAX = CABS1( AP( KX ) )
                     JMAX = J
                  }
                  KX = KX + J
               } // 20
               KPC = ( IMAX-1 )*IMAX / 2 + 1
               if ( IMAX.GT.1 ) {
                  JMAX = IZAMAX( IMAX-1, AP( KPC ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( AP( KPC+JMAX-1 ) ) )
               }

               if ( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K
               } else if ( ABS( DBLE( AP( KPC+IMAX-1 ) ) ).GE.ALPHA* ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX
               } else {

                  // interchange rows and columns K-1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX
                  KSTEP = 2
               }
            }

            KK = K - KSTEP + 1
            IF( KSTEP.EQ.2 ) KNC = KNC - K + 1
            if ( KP.NE.KK ) {

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               zswap(KP-1, AP( KNC ), 1, AP( KPC ), 1 );
               KX = KPC + KP - 1
               DO 30 J = KP + 1, KK - 1
                  KX = KX + J - 1
                  T = DCONJG( AP( KNC+J-1 ) )
                  AP( KNC+J-1 ) = DCONJG( AP( KX ) )
                  AP( KX ) = T
               } // 30
               AP( KX+KK-1 ) = DCONJG( AP( KX+KK-1 ) )
               R1 = DBLE( AP( KNC+KK-1 ) )
               AP( KNC+KK-1 ) = DBLE( AP( KPC+KP-1 ) )
               AP( KPC+KP-1 ) = R1
               if ( KSTEP.EQ.2 ) {
                  AP( KC+K-1 ) = DBLE( AP( KC+K-1 ) )
                  T = AP( KC+K-2 )
                  AP( KC+K-2 ) = AP( KC+KP-1 )
                  AP( KC+KP-1 ) = T
               }
            } else {
               AP( KC+K-1 ) = DBLE( AP( KC+K-1 ) )
               IF( KSTEP.EQ.2 ) AP( KC-1 ) = DBLE( AP( KC-1 ) )
            }

            // Update the leading submatrix

            if ( KSTEP.EQ.1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Perform a rank-1 update of A(1:k-1,1:k-1) as

               // A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H

               R1 = ONE / DBLE( AP( KC+K-1 ) )
               zhpr(UPLO, K-1, -R1, AP( KC ), 1, AP );

               // Store U(k) in column k

               zdscal(K-1, R1, AP( KC ), 1 );
            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
                  // = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H

               if ( K.GT.2 ) {

                  D = DLAPY2( DBLE( AP( K-1+( K-1 )*K / 2 ) ), DIMAG( AP( K-1+( K-1 )*K / 2 ) ) )
                  D22 = DBLE( AP( K-1+( K-2 )*( K-1 ) / 2 ) ) / D
                  D11 = DBLE( AP( K+( K-1 )*K / 2 ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D12 = AP( K-1+( K-1 )*K / 2 ) / D
                  D = TT / D

                  DO 50 J = K - 2, 1, -1
                     WKM1 = D*( D11*AP( J+( K-2 )*( K-1 ) / 2 )- DCONJG( D12 )*AP( J+( K-1 )*K / 2 ) )                      WK = D*( D22*AP( J+( K-1 )*K / 2 )-D12* AP( J+( K-2 )*( K-1 ) / 2 ) )
                     DO 40 I = J, 1, -1
                        AP( I+( J-1 )*J / 2 ) = AP( I+( J-1 )*J / 2 ) - AP( I+( K-1 )*K / 2 )*DCONJG( WK ) - AP( I+( K-2 )*( K-1 ) / 2 )*DCONJG( WKM1 )
                     } // 40
                     AP( J+( K-1 )*K / 2 ) = WK
                     AP( J+( K-2 )*( K-1 ) / 2 ) = WKM1
                     AP( J+( J-1 )*J / 2 ) = DCMPLX( DBLE( AP( J+( J- 1 )*J / 2 ) ), 0.0D+0 )
                  } // 50

               }

            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP.EQ.1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP
         KC = KNC - K
         GO TO 10

      } else {

         // Factorize A as L*D*L**H using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1
         KC = 1
         NPP = N*( N+1 ) / 2
         } // 60
         KNC = KC

         // If K > N, exit from loop

         IF( K.GT.N ) GO TO 110
         KSTEP = 1

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( DBLE( AP( KC ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value

         if ( K.LT.N ) {
            IMAX = K + IZAMAX( N-K, AP( KC+1 ), 1 )
            COLMAX = CABS1( AP( KC+IMAX-K ) )
         } else {
            COLMAX = ZERO
         }

         if ( MAX( ABSAKK, COLMAX ).EQ.ZERO ) {

            // Column K is zero: set INFO and continue

            IF( INFO.EQ.0 ) INFO = K
            KP = K
            AP( KC ) = DBLE( AP( KC ) )
         } else {
            if ( ABSAKK.GE.ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               ROWMAX = ZERO
               KX = KC + IMAX - K
               DO 70 J = K, IMAX - 1
                  if ( CABS1( AP( KX ) ).GT.ROWMAX ) {
                     ROWMAX = CABS1( AP( KX ) )
                     JMAX = J
                  }
                  KX = KX + N - J
               } // 70
               KPC = NPP - ( N-IMAX+1 )*( N-IMAX+2 ) / 2 + 1
               if ( IMAX.LT.N ) {
                  JMAX = IMAX + IZAMAX( N-IMAX, AP( KPC+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( AP( KPC+JMAX-IMAX ) ) )
               }

               if ( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K
               } else if ( ABS( DBLE( AP( KPC ) ) ).GE.ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX
               } else {

                  // interchange rows and columns K+1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX
                  KSTEP = 2
               }
            }

            KK = K + KSTEP - 1
            IF( KSTEP.EQ.2 ) KNC = KNC + N - K + 1
            if ( KP.NE.KK ) {

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               IF( KP.LT.N ) CALL ZSWAP( N-KP, AP( KNC+KP-KK+1 ), 1, AP( KPC+1 ), 1 )
               KX = KNC + KP - KK
               DO 80 J = KK + 1, KP - 1
                  KX = KX + N - J + 1
                  T = DCONJG( AP( KNC+J-KK ) )
                  AP( KNC+J-KK ) = DCONJG( AP( KX ) )
                  AP( KX ) = T
               } // 80
               AP( KNC+KP-KK ) = DCONJG( AP( KNC+KP-KK ) )
               R1 = DBLE( AP( KNC ) )
               AP( KNC ) = DBLE( AP( KPC ) )
               AP( KPC ) = R1
               if ( KSTEP.EQ.2 ) {
                  AP( KC ) = DBLE( AP( KC ) )
                  T = AP( KC+1 )
                  AP( KC+1 ) = AP( KC+KP-K )
                  AP( KC+KP-K ) = T
               }
            } else {
               AP( KC ) = DBLE( AP( KC ) )
               IF( KSTEP.EQ.2 ) AP( KNC ) = DBLE( AP( KNC ) )
            }

            // Update the trailing submatrix

            if ( KSTEP.EQ.1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               if ( K.LT.N ) {

                  // Perform a rank-1 update of A(k+1:n,k+1:n) as

                  // A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H

                  R1 = ONE / DBLE( AP( KC ) )
                  zhpr(UPLO, N-K, -R1, AP( KC+1 ), 1, AP( KC+N-K+1 ) );

                  // Store L(k) in column K

                  zdscal(N-K, R1, AP( KC+1 ), 1 );
               }
            } else {

               // 2-by-2 pivot block D(k): columns K and K+1 now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               if ( K.LT.N-1 ) {

                  // Perform a rank-2 update of A(k+2:n,k+2:n) as

                  // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
                     // = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H

                  // where L(k) and L(k+1) are the k-th and (k+1)-th
                  // columns of L

                  D = DLAPY2( DBLE( AP( K+1+( K-1 )*( 2*N-K ) / 2 ) ), DIMAG( AP( K+1+( K-1 )*( 2*N-K ) / 2 ) ) )
                  D11 = DBLE( AP( K+1+K*( 2*N-K-1 ) / 2 ) ) / D
                  D22 = DBLE( AP( K+( K-1 )*( 2*N-K ) / 2 ) ) / D
                  TT = ONE / ( D11*D22-ONE )
                  D21 = AP( K+1+( K-1 )*( 2*N-K ) / 2 ) / D
                  D = TT / D

                  DO 100 J = K + 2, N
                     WK = D*( D11*AP( J+( K-1 )*( 2*N-K ) / 2 )-D21* AP( J+K*( 2*N-K-1 ) / 2 ) )                      WKP1 = D*( D22*AP( J+K*( 2*N-K-1 ) / 2 )- DCONJG( D21 )*AP( J+( K-1 )*( 2*N-K ) / 2 ) )
                     for (I = J; I <= N; I++) { // 90
                        AP( I+( J-1 )*( 2*N-J ) / 2 ) = AP( I+( J-1 )* ( 2*N-J ) / 2 ) - AP( I+( K-1 )*( 2*N-K ) / 2 )*DCONJG( WK ) - AP( I+K*( 2*N-K-1 ) / 2 )* DCONJG( WKP1 )
                     } // 90
                     AP( J+( K-1 )*( 2*N-K ) / 2 ) = WK
                     AP( J+K*( 2*N-K-1 ) / 2 ) = WKP1
                     AP( J+( J-1 )*( 2*N-J ) / 2 ) = DCMPLX( DBLE( AP( J+( J-1 )*( 2*N-J ) / 2 ) ), 0.0D+0 )
                  } // 100
               }
            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP.EQ.1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP
         KC = KNC + N - K + 2
         GO TO 60

      }

      } // 110
      RETURN

      // End of ZHPTRF

      }

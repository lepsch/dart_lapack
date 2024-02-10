      void csptrf(UPLO, N, AP, IPIV, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N;
      int                IPIV( * );
      Complex            AP( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      int                I, IMAX, J, JMAX, K, KC, KK, KNC, KP, KPC, KSTEP, KX, NPP;
      double               ABSAKK, ALPHA, COLMAX, ROWMAX;
      Complex            D11, D12, D21, D22, R1, T, WK, WKM1, WKP1, ZDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ICAMAX;
      // EXTERNAL lsame, ICAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSPR, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('CSPTRF', -INFO );
         return;
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      if ( UPPER ) {

         // Factorize A as U*D*U**T using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N;
         KC = ( N-1 )*N / 2 + 1;
         } // 10
         KNC = KC;

         // If K < 1, exit from loop

         if (K < 1) GO TO 110;
         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = CABS1( AP( KC+K-1 ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value

         if ( K > 1 ) {
            IMAX = ICAMAX( K-1, AP( KC ), 1 );
            COLMAX = CABS1( AP( KC+IMAX-1 ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
         } else {
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               ROWMAX = ZERO;
               JMAX = IMAX;
               KX = IMAX*( IMAX+1 ) / 2 + IMAX;
               for (J = IMAX + 1; J <= K; J++) { // 20
                  if ( CABS1( AP( KX ) ) > ROWMAX ) {
                     ROWMAX = CABS1( AP( KX ) );
                     JMAX = J;
                  }
                  KX = KX + J;
               } // 20
               KPC = ( IMAX-1 )*IMAX / 2 + 1;
               if ( IMAX > 1 ) {
                  JMAX = ICAMAX( IMAX-1, AP( KPC ), 1 );
                  ROWMAX = max( ROWMAX, CABS1( AP( KPC+JMAX-1 ) ) );
               }

               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;
               } else if ( CABS1( AP( KPC+IMAX-1 ) ) >= ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX;
               } else {

                  // interchange rows and columns K-1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX;
                  KSTEP = 2;
               }
            }

            KK = K - KSTEP + 1;
            if (KSTEP == 2) KNC = KNC - K + 1;
            if ( KP != KK ) {

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               cswap(KP-1, AP( KNC ), 1, AP( KPC ), 1 );
               KX = KPC + KP - 1;
               for (J = KP + 1; J <= KK - 1; J++) { // 30
                  KX = KX + J - 1;
                  T = AP( KNC+J-1 );
                  AP[KNC+J-1] = AP( KX );
                  AP[KX] = T;
               } // 30
               T = AP( KNC+KK-1 );
               AP[KNC+KK-1] = AP( KPC+KP-1 );
               AP[KPC+KP-1] = T;
               if ( KSTEP == 2 ) {
                  T = AP( KC+K-2 );
                  AP[KC+K-2] = AP( KC+KP-1 );
                  AP[KC+KP-1] = T;
               }
            }

            // Update the leading submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Perform a rank-1 update of A(1:k-1,1:k-1) as

               // A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T

               R1 = CONE / AP( KC+K-1 );
               cspr(UPLO, K-1, -R1, AP( KC ), 1, AP );

               // Store U(k) in column k

               cscal(K-1, R1, AP( KC ), 1 );
            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
                  // = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T

               if ( K > 2 ) {

                  D12 = AP( K-1+( K-1 )*K / 2 );
                  D22 = AP( K-1+( K-2 )*( K-1 ) / 2 ) / D12;
                  D11 = AP( K+( K-1 )*K / 2 ) / D12;
                  T = CONE / ( D11*D22-CONE );
                  D12 = T / D12;

                  for (J = K - 2; J >= 1; J--) { // 50
                     WKM1 = D12*( D11*AP( J+( K-2 )*( K-1 ) / 2 )- AP( J+( K-1 )*K / 2 ) )                      WK = D12*( D22*AP( J+( K-1 )*K / 2 )- AP( J+( K-2 )*( K-1 ) / 2 ) );
                     for (I = J; I >= 1; I--) { // 40
                        AP[I+( J-1 )*J / 2] = AP( I+( J-1 )*J / 2 ) - AP( I+( K-1 )*K / 2 )*WK - AP( I+( K-2 )*( K-1 ) / 2 )*WKM1;
                     } // 40
                     AP[J+( K-1 )*K / 2] = WK;
                     AP[J+( K-2 )*( K-1 ) / 2] = WKM1;
                  } // 50

               }
            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -KP;
            IPIV[K-1] = -KP;
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP;
         KC = KNC - K;
         GO TO 10;

      } else {

         // Factorize A as L*D*L**T using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1;
         KC = 1;
         NPP = N*( N+1 ) / 2;
         } // 60
         KNC = KC;

         // If K > N, exit from loop

         if (K > N) GO TO 110;
         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = CABS1( AP( KC ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value

         if ( K < N ) {
            IMAX = K + ICAMAX( N-K, AP( KC+1 ), 1 );
            COLMAX = CABS1( AP( KC+IMAX-K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
         } else {
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               ROWMAX = ZERO;
               KX = KC + IMAX - K;
               for (J = K; J <= IMAX - 1; J++) { // 70
                  if ( CABS1( AP( KX ) ) > ROWMAX ) {
                     ROWMAX = CABS1( AP( KX ) );
                     JMAX = J;
                  }
                  KX = KX + N - J;
               } // 70
               KPC = NPP - ( N-IMAX+1 )*( N-IMAX+2 ) / 2 + 1;
               if ( IMAX < N ) {
                  JMAX = IMAX + ICAMAX( N-IMAX, AP( KPC+1 ), 1 );
                  ROWMAX = max( ROWMAX, CABS1( AP( KPC+JMAX-IMAX ) ) );
               }

               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;
               } else if ( CABS1( AP( KPC ) ) >= ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX;
               } else {

                  // interchange rows and columns K+1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX;
                  KSTEP = 2;
               }
            }

            KK = K + KSTEP - 1;
            if (KSTEP == 2) KNC = KNC + N - K + 1;
            if ( KP != KK ) {

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               if (KP < N) cswap( N-KP, AP( KNC+KP-KK+1 ), 1, AP( KPC+1 ), 1 );
               KX = KNC + KP - KK;
               for (J = KK + 1; J <= KP - 1; J++) { // 80
                  KX = KX + N - J + 1;
                  T = AP( KNC+J-KK );
                  AP[KNC+J-KK] = AP( KX );
                  AP[KX] = T;
               } // 80
               T = AP( KNC );
               AP[KNC] = AP( KPC );
               AP[KPC] = T;
               if ( KSTEP == 2 ) {
                  T = AP( KC+1 );
                  AP[KC+1] = AP( KC+KP-K );
                  AP[KC+KP-K] = T;
               }
            }

            // Update the trailing submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               if ( K < N ) {

                  // Perform a rank-1 update of A(k+1:n,k+1:n) as

                  // A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T

                  R1 = CONE / AP( KC );
                  cspr(UPLO, N-K, -R1, AP( KC+1 ), 1, AP( KC+N-K+1 ) );

                  // Store L(k) in column K

                  cscal(N-K, R1, AP( KC+1 ), 1 );
               }
            } else {

               // 2-by-2 pivot block D(k): columns K and K+1 now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               if ( K < N-1 ) {

                  // Perform a rank-2 update of A(k+2:n,k+2:n) as

                  // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T
                     // = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T

                  // where L(k) and L(k+1) are the k-th and (k+1)-th
                  // columns of L

                  D21 = AP( K+1+( K-1 )*( 2*N-K ) / 2 );
                  D11 = AP( K+1+K*( 2*N-K-1 ) / 2 ) / D21;
                  D22 = AP( K+( K-1 )*( 2*N-K ) / 2 ) / D21;
                  T = CONE / ( D11*D22-CONE );
                  D21 = T / D21;

                  for (J = K + 2; J <= N; J++) { // 100
                     WK = D21*( D11*AP( J+( K-1 )*( 2*N-K ) / 2 )- AP( J+K*( 2*N-K-1 ) / 2 ) )                      WKP1 = D21*( D22*AP( J+K*( 2*N-K-1 ) / 2 )- AP( J+( K-1 )*( 2*N-K ) / 2 ) );
                     for (I = J; I <= N; I++) { // 90
                        AP[I+( J-1 )*( 2*N-J ) / 2] = AP( I+( J-1 )* ( 2*N-J ) / 2 ) - AP( I+( K-1 )*( 2*N-K ) / 2 )*WK - AP( I+K*( 2*N-K-1 ) / 2 )*WKP1;
                     } // 90
                     AP[J+( K-1 )*( 2*N-K ) / 2] = WK;
                     AP[J+K*( 2*N-K-1 ) / 2] = WKP1;
                  } // 100
               }
            }
         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -KP;
            IPIV[K+1] = -KP;
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP;
         KC = KNC + N - K + 2;
         GO TO 60;

      }

      } // 110
      }

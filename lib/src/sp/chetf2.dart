      void chetf2(UPLO, N, A, LDA, IPIV, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IMAX, J, JMAX, K, KK, KP, KSTEP;
      double               ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, ROWMAX, TT;
      Complex            D12, D21, T, WK, WKM1, WKP1, ZDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      //- int                ICAMAX;
      //- REAL               SLAPY2;
      // EXTERNAL lsame, ICAMAX, SLAPY2, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CSSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CHETF2', -INFO );
         return;
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      if ( UPPER ) {

         // Factorize A as U*D*U**H using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N;
         } // 10

         // If K < 1, exit from loop

         if (K < 1) GO TO 90;
         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( double( A( K, K ) ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = ICAMAX( K-1, A( 1, K ), 1 );
            COLMAX = CABS1( A( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( (max( ABSAKK, COLMAX ) == ZERO) || SISNAN(ABSAKK) ) {

            // Column K is or underflow, or contains a NaN:
            // set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            A[K][K] = double( A( K, K ) );
         } else {
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               JMAX = IMAX + ICAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA );
               ROWMAX = CABS1( A( IMAX, JMAX ) );
               if ( IMAX > 1 ) {
                  JMAX = ICAMAX( IMAX-1, A( 1, IMAX ), 1 );
                  ROWMAX = max( ROWMAX, CABS1( A( JMAX, IMAX ) ) );
               }

               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;
               } else if ( ABS( double( A( IMAX, IMAX ) ) ) >= ALPHA*ROWMAX ) {

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
            if ( KP != KK ) {

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               cswap(KP-1, A( 1, KK ), 1, A( 1, KP ), 1 );
               for (J = KP + 1; J <= KK - 1; J++) { // 20
                  T = CONJG( A( J, KK ) );
                  A[J][KK] = CONJG( A( KP, J ) );
                  A[KP][J] = T;
               } // 20
               A[KP][KK] = CONJG( A( KP, KK ) );
               R1 = double( A( KK, KK ) );
               A[KK][KK] = double( A( KP, KP ) );
               A[KP][KP] = R1;
               if ( KSTEP == 2 ) {
                  A[K][K] = double( A( K, K ) );
                  T = A( K-1, K );
                  A[K-1, K] = A( KP, K );
                  A[KP][K] = T;
               }
            } else {
               A[K][K] = double( A( K, K ) );
               if (KSTEP == 2) A( K-1, K-1 ) = double( A( K-1, K-1 ) );
            }

            // Update the leading submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Perform a rank-1 update of A(1:k-1,1:k-1) as

               // A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H

               R1 = ONE / REAL( A( K, K ) );
               cher(UPLO, K-1, -R1, A( 1, K ), 1, A, LDA );

               // Store U(k) in column k

               csscal(K-1, R1, A( 1, K ), 1 );
            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
                  // = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H

               if ( K > 2 ) {

                  D = SLAPY2( double( A( K-1, K ) ), AIMAG( A( K-1, K ) ) );
                  D22 = double( A( K-1, K-1 ) ) / D;
                  D11 = double( A( K, K ) ) / D;
                  TT = ONE / ( D11*D22-ONE );
                  D12 = A( K-1, K ) / D;
                  D = TT / D;

                  for (J = K - 2; J >= 1; J--) { // 40
                     WKM1 = D*( D11*A( J, K-1 )-CONJG( D12 )*A( J, K ) );
                     WK = D*( D22*A( J, K )-D12*A( J, K-1 ) );
                     for (I = J; I >= 1; I--) { // 30
                        A[I][J] = A( I, J ) - A( I, K )*CONJG( WK ) - A( I, K-1 )*CONJG( WKM1 );
                     } // 30
                     A[J][K] = WK;
                     A[J, K-1] = WKM1;
                     A[J][J] = CMPLX( double( A( J, J ) ), 0.0 );
                  } // 40

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
         GO TO 10;

      } else {

         // Factorize A as L*D*L**H using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1;
         } // 50

         // If K > N, exit from loop

         if (K > N) GO TO 90;
         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( double( A( K, K ) ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + ICAMAX( N-K, A( K+1, K ), 1 );
            COLMAX = CABS1( A( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( (max( ABSAKK, COLMAX ) == ZERO) || SISNAN(ABSAKK) ) {

            // Column K is zero or underflow, contains a NaN:
            // set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            A[K][K] = double( A( K, K ) );
         } else {
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               JMAX = K - 1 + ICAMAX( IMAX-K, A( IMAX, K ), LDA );
               ROWMAX = CABS1( A( IMAX, JMAX ) );
               if ( IMAX < N ) {
                  JMAX = IMAX + ICAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 );
                  ROWMAX = max( ROWMAX, CABS1( A( JMAX, IMAX ) ) );
               }

               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;
               } else if ( ABS( double( A( IMAX, IMAX ) ) ) >= ALPHA*ROWMAX ) {

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
            if ( KP != KK ) {

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               if (KP < N) cswap( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 );
               for (J = KK + 1; J <= KP - 1; J++) { // 60
                  T = CONJG( A( J, KK ) );
                  A[J][KK] = CONJG( A( KP, J ) );
                  A[KP][J] = T;
               } // 60
               A[KP][KK] = CONJG( A( KP, KK ) );
               R1 = double( A( KK, KK ) );
               A[KK][KK] = double( A( KP, KP ) );
               A[KP][KP] = R1;
               if ( KSTEP == 2 ) {
                  A[K][K] = double( A( K, K ) );
                  T = A( K+1, K );
                  A[K+1, K] = A( KP, K );
                  A[KP][K] = T;
               }
            } else {
               A[K][K] = double( A( K, K ) );
               if (KSTEP == 2) A( K+1, K+1 ) = double( A( K+1, K+1 ) );
            }

            // Update the trailing submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               if ( K < N ) {

                  // Perform a rank-1 update of A(k+1:n,k+1:n) as

                  // A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H

                  R1 = ONE / REAL( A( K, K ) );
                  cher(UPLO, N-K, -R1, A( K+1, K ), 1, A( K+1, K+1 ), LDA );

                  // Store L(k) in column K

                  csscal(N-K, R1, A( K+1, K ), 1 );
               }
            } else {

               // 2-by-2 pivot block D(k)

               if ( K < N-1 ) {

                  // Perform a rank-2 update of A(k+2:n,k+2:n) as

                  // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
                     // = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H

                  // where L(k) and L(k+1) are the k-th and (k+1)-th
                  // columns of L

                  D = SLAPY2( double( A( K+1, K ) ), AIMAG( A( K+1, K ) ) );
                  D11 = double( A( K+1, K+1 ) ) / D;
                  D22 = double( A( K, K ) ) / D;
                  TT = ONE / ( D11*D22-ONE );
                  D21 = A( K+1, K ) / D;
                  D =  TT / D;

                  for (J = K + 2; J <= N; J++) { // 80
                     WK = D*( D11*A( J, K )-D21*A( J, K+1 ) );
                     WKP1 = D*( D22*A( J, K+1 )-CONJG( D21 )*A( J, K ) );
                     for (I = J; I <= N; I++) { // 70
                        A[I][J] = A( I, J ) - A( I, K )*CONJG( WK ) - A( I, K+1 )*CONJG( WKP1 );
                     } // 70
                     A[J][K] = WK;
                     A[J, K+1] = WKP1;
                     A[J][J] = CMPLX( double( A( J, J ) ), 0.0 );
                  } // 80
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
         GO TO 50;

      }

      } // 90
      return;
      }

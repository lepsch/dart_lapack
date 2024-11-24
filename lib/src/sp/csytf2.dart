// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void csytf2(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      int                IPIV( * );
      Complex            A( LDA, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               UPPER;
      int                I, IMAX, J, JMAX, K, KK, KP, KSTEP;
      double               ABSAKK, ALPHA, COLMAX, ROWMAX;
      Complex            D11, D12, D21, D22, R1, T, WK, WKM1, WKP1, Z;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      //- int                ICAMAX;
      // EXTERNAL lsame, ICAMAX, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSCAL, CSWAP, CSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[Z] = ( double( Z ) ).abs() + ( AIMAG( Z ) ).abs();

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
         xerbla('CSYTF2', -INFO );
         return;
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      if ( UPPER ) {

         // Factorize A as U*D*U**T using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N;
         } // 10

         // If K < 1, exit from loop

         if (K < 1) GO TO 70;
         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = CABS1( A( K, K ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = ICAMAX( K-1, A( 1, K ), 1 );
            COLMAX = CABS1( A( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO || SISNAN(ABSAKK) ) {

            // Column K is zero or underflow, or contains a NaN:
            // set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
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
               } else if ( CABS1( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX ) {

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
               cswap(KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA );
               T = A( KK, KK );
               A[KK][KK] = A( KP, KP );
               A[KP][KP] = T;
               if ( KSTEP == 2 ) {
                  T = A( K-1, K );
                  A[K-1][K] = A( KP, K );
                  A[KP][K] = T;
               }
            }

            // Update the leading submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               // Perform a rank-1 update of A(1:k-1,1:k-1) as

               // A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T

               R1 = CONE / A( K, K );
               csyr(UPLO, K-1, -R1, A( 1, K ), 1, A, LDA );

               // Store U(k) in column k

               cscal(K-1, R1, A( 1, K ), 1 );
            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
               //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T

               if ( K > 2 ) {

                  D12 = A( K-1, K );
                  D22 = A( K-1, K-1 ) / D12;
                  D11 = A( K, K ) / D12;
                  T = CONE / ( D11*D22-CONE );
                  D12 = T / D12;

                  for (J = K - 2; J >= 1; J--) { // 30
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) );
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) );
                     for (I = J; I >= 1; I--) { // 20
                        A[I][J] = A( I, J ) - A( I, K )*WK - A( I, K-1 )*WKM1;
                     } // 20
                     A[J][K] = WK;
                     A[J][K-1] = WKM1;
                  } // 30

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

         // Factorize A as L*D*L**T using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1;
         } // 40

         // If K > N, exit from loop

         if (K > N) GO TO 70;
         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = CABS1( A( K, K ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + ICAMAX( N-K, A( K+1, K ), 1 );
            COLMAX = CABS1( A( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO || SISNAN(ABSAKK) ) {

            // Column K is zero or underflow, or contains a NaN:
            // set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
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
               } else if ( CABS1( A( IMAX, IMAX ) ) >= ALPHA*ROWMAX ) {

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
               cswap(KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA );
               T = A( KK, KK );
               A[KK][KK] = A( KP, KP );
               A[KP][KP] = T;
               if ( KSTEP == 2 ) {
                  T = A( K+1, K );
                  A[K+1][K] = A( KP, K );
                  A[KP][K] = T;
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

                  R1 = CONE / A( K, K );
                  csyr(UPLO, N-K, -R1, A( K+1, K ), 1, A( K+1, K+1 ), LDA );

                  // Store L(k) in column K

                  cscal(N-K, R1, A( K+1, K ), 1 );
               }
            } else {

               // 2-by-2 pivot block D(k)

               if ( K < N-1 ) {

                  // Perform a rank-2 update of A(k+2:n,k+2:n) as

                  // A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T
                  //    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T

                  // where L(k) and L(k+1) are the k-th and (k+1)-th
                  // columns of L

                  D21 = A( K+1, K );
                  D11 = A( K+1, K+1 ) / D21;
                  D22 = A( K, K ) / D21;
                  T = CONE / ( D11*D22-CONE );
                  D21 = T / D21;

                  for (J = K + 2; J <= N; J++) { // 60
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) );
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) );
                     for (I = J; I <= N; I++) { // 50
                        A[I][J] = A( I, J ) - A( I, K )*WK - A( I, K+1 )*WKP1;
                     } // 50
                     A[J][K] = WK;
                     A[J][K+1] = WKP1;
                  } // 60
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
         GO TO 40;

      }

      } // 70
      }

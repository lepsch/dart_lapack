// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void chetri_rook(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      int                IPIV( * );
      Complex            A( LDA, * ), WORK( * );
      // ..

      double               ONE;
      Complex            CONE, CZERO;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      bool               UPPER;
      int                J, K, KP, KSTEP;
      double               AK, AKP1, D, T;
      Complex            AKKP1, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHEMV, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL

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
         xerbla('CHETRI_ROOK', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         for (INFO = N; INFO >= 1; INFO--) { // 10
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == CZERO ) return;
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) { // 20
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == CZERO ) return;
         } // 20
      }
      INFO = 0;

      if ( UPPER ) {

         // Compute inv(A) from the factorization A = U*D*U**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1;
         } // 30

         // If K > N, exit from loop.

         if (K > N) GO TO 70;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A[K][K] = ONE / REAL( A( K, K ) );

            // Compute column K of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               chemv[UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K][K] = A( K, K ) - double( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ( A( K, K+1 ) ).abs();
            AK = double( A( K, K ) ) / T;
            AKP1 = double( A( K+1, K+1 ) ) / T;
            AKKP1 = A( K, K+1 ) / T;
            D = T*( AK*AKP1-ONE );
            A[K][K] = AKP1 / D;
            A[K+1][K+1] = AK / D;
            A[K][K+1] = -AKKP1 / D;

            // Compute columns K and K+1 of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               chemv[UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K][K] = A( K, K ) - double( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )                A( K, K+1 ) = A( K, K+1 ) - CDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 );
               ccopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               chemv[UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K+1 ), 1 )                A( K+1, K+1] = A( K+1, K+1 ) - double( CDOTC( K-1, WORK, 1, A( 1, K+1 ), 1 ) );
            }
            KSTEP = 2;
         }

         if ( KSTEP == 1 ) {

            // Interchange rows and columns K and IPIV(K) in the leading
            // submatrix A(1:k,1:k)

            KP = IPIV( K );
            if ( KP != K ) {

               if (KP > 1) cswap( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );

               for (J = KP + 1; J <= K - 1; J++) { // 40
                  TEMP = CONJG( A( J, K ) );
                  A[J][K] = CONJG( A( KP, J ) );
                  A[KP][J] = TEMP;
               } // 40

               A[KP][K] = CONJG( A( KP, K ) );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         } else {

            // Interchange rows and columns K and K+1 with -IPIV(K) and
            // -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n)

            // (1) Interchange rows and columns K and -IPIV(K)

            KP = -IPIV( K );
            if ( KP != K ) {

               if (KP > 1) cswap( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );

               for (J = KP + 1; J <= K - 1; J++) { // 50
                  TEMP = CONJG( A( J, K ) );
                  A[J][K] = CONJG( A( KP, J ) );
                  A[KP][J] = TEMP;
               } // 50

               A[KP][K] = CONJG( A( KP, K ) );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;

               TEMP = A( K, K+1 );
               A[K][K+1] = A( KP, K+1 );
               A[KP][K+1] = TEMP;
            }

            // (2) Interchange rows and columns K+1 and -IPIV(K+1)

            K = K + 1;
            KP = -IPIV( K );
            if ( KP != K ) {

               if (KP > 1) cswap( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );

               for (J = KP + 1; J <= K - 1; J++) { // 60
                  TEMP = CONJG( A( J, K ) );
                  A[J][K] = CONJG( A( KP, J ) );
                  A[KP][J] = TEMP;
               } // 60

               A[KP][K] = CONJG( A( KP, K ) );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         }

         K = K + 1;
         GO TO 30;
         } // 70

      } else {

         // Compute inv(A) from the factorization A = L*D*L**H.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N;
         } // 80

         // If K < 1, exit from loop.

         if (K < 1) GO TO 120;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A[K][K] = ONE / REAL( A( K, K ) );

            // Compute column K of the inverse.

            if ( K < N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               chemv[UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K][K] = A( K, K ) - double( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ( A( K, K-1 ) ).abs();
            AK = double( A( K-1, K-1 ) ) / T;
            AKP1 = double( A( K, K ) ) / T;
            AKKP1 = A( K, K-1 ) / T;
            D = T*( AK*AKP1-ONE );
            A[K-1][K-1] = AKP1 / D;
            A[K][K] = AK / D;
            A[K][K-1] = -AKKP1 / D;

            // Compute columns K-1 and K of the inverse.

            if ( K < N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               chemv[UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K][K] = A( K, K ) - double( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )                A( K, K-1 ) = A( K, K-1 ) - CDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               ccopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               chemv[UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1] = A( K-1, K-1 ) - double( CDOTC( N-K, WORK, 1, A( K+1, K-1 ), 1 ) );
            }
            KSTEP = 2;
         }

         if ( KSTEP == 1 ) {

            // Interchange rows and columns K and IPIV(K) in the trailing
            // submatrix A(k:n,k:n)

            KP = IPIV( K );
            if ( KP != K ) {

               if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );

               for (J = K + 1; J <= KP - 1; J++) { // 90
                  TEMP = CONJG( A( J, K ) );
                  A[J][K] = CONJG( A( KP, J ) );
                  A[KP][J] = TEMP;
               } // 90

               A[KP][K] = CONJG( A( KP, K ) );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         } else {

            // Interchange rows and columns K and K-1 with -IPIV(K) and
            // -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)

            // (1) Interchange rows and columns K and -IPIV(K)

            KP = -IPIV( K );
            if ( KP != K ) {

               if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );

               for (J = K + 1; J <= KP - 1; J++) { // 100
                  TEMP = CONJG( A( J, K ) );
                  A[J][K] = CONJG( A( KP, J ) );
                  A[KP][J] = TEMP;
              } // 100

               A[KP][K] = CONJG( A( KP, K ) );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;

               TEMP = A( K, K-1 );
               A[K][K-1] = A( KP, K-1 );
               A[KP][K-1] = TEMP;
            }

            // (2) Interchange rows and columns K-1 and -IPIV(K-1)

            K = K - 1;
            KP = -IPIV( K );
            if ( KP != K ) {

               if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );

               for (J = K + 1; J <= KP - 1; J++) { // 110
                  TEMP = CONJG( A( J, K ) );
                  A[J][K] = CONJG( A( KP, J ) );
                  A[KP][J] = TEMP;
              } // 110

               A[KP][K] = CONJG( A( KP, K ) );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         }

         K = K - 1;
         GO TO 80;
         } // 120
      }

      }

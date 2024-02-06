      void csytri_rook(UPLO, N, A, LDA, IPIV, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex            A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                K, KP, KSTEP;
      Complex            AK, AKKP1, AKP1, D, T, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- COMPLEX            CDOTU;
      // EXTERNAL lsame, CDOTU
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CSWAP, CSYMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
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
         xerbla('CSYTRI_ROOK', -INFO );
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

         // Compute inv(A) from the factorization A = U*D*U**T.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1;
         } // 30

         // If K > N, exit from loop.

         if (K > N) GO TO 40;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A[K][K] = CONE / A( K, K );

            // Compute column K of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               csymv[UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K, K] = A( K, K ) - CDOTU( K-1, WORK, 1, A( 1, K ), 1 );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = A( K, K+1 );
            AK = A( K, K ) / T;
            AKP1 = A( K+1, K+1 ) / T;
            AKKP1 = A( K, K+1 ) / T;
            D = T*( AK*AKP1-CONE );
            A[K][K] = AKP1 / D;
            A[K+1, K+1] = AK / D;
            A[K, K+1] = -AKKP1 / D;

            // Compute columns K and K+1 of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               csymv[UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K, K] = A( K, K ) - CDOTU( K-1, WORK, 1, A( 1, K ), 1 )                A( K, K+1 ) = A( K, K+1 ) - CDOTU( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 );
               ccopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               csymv[UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K+1 ), 1 )                A( K+1, K+1] = A( K+1, K+1 ) - CDOTU( K-1, WORK, 1, A( 1, K+1 ), 1 );
            }
            KSTEP = 2;
         }

         if ( KSTEP == 1 ) {

            // Interchange rows and columns K and IPIV(K) in the leading
            // submatrix A(1:k+1,1:k+1)

            KP = IPIV( K );
            if ( KP != K ) {
               if (KP > 1) cswap( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
               cswap(K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA );
               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         } else {

            // Interchange rows and columns K and K+1 with -IPIV(K) and
            // -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1)

            KP = -IPIV( K );
            if ( KP != K ) {
               if (KP > 1) cswap( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
               cswap(K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
               TEMP = A( K, K+1 );
               A[K, K+1] = A( KP, K+1 );
               A[KP, K+1] = TEMP;
            }

            K = K + 1;
            KP = -IPIV( K );
            if ( KP != K ) {
               if (KP > 1) cswap( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
               cswap(K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA );
               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         }

         K = K + 1;
         GO TO 30;
         } // 40

      } else {

         // Compute inv(A) from the factorization A = L*D*L**T.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N;
         } // 50

         // If K < 1, exit from loop.

         if (K < 1) GO TO 60;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A[K][K] = CONE / A( K, K );

            // Compute column K of the inverse.

            if ( K < N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               csymv[UPLO, N-K,-CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K, K] = A( K, K ) - CDOTU( N-K, WORK, 1, A( K+1, K ), 1 );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = A( K, K-1 );
            AK = A( K-1, K-1 ) / T;
            AKP1 = A( K, K ) / T;
            AKKP1 = A( K, K-1 ) / T;
            D = T*( AK*AKP1-CONE );
            A[K-1, K-1] = AKP1 / D;
            A[K][K] = AK / D;
            A[K, K-1] = -AKKP1 / D;

            // Compute columns K-1 and K of the inverse.

            if ( K < N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               csymv[UPLO, N-K,-CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K, K] = A( K, K ) - CDOTU( N-K, WORK, 1, A( K+1, K ), 1 )                A( K, K-1 ) = A( K, K-1 ) - CDOTU( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               ccopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               csymv[UPLO, N-K,-CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1] = A( K-1, K-1 ) - CDOTU( N-K, WORK, 1, A( K+1, K-1 ), 1 );
            }
            KSTEP = 2;
         }

         if ( KSTEP == 1 ) {

            // Interchange rows and columns K and IPIV(K) in the trailing
            // submatrix A(k-1:n,k-1:n)

            KP = IPIV( K );
            if ( KP != K ) {
               if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );
               cswap(KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA );
               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         } else {

            // Interchange rows and columns K and K-1 with -IPIV(K) and
            // -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)

            KP = -IPIV( K );
            if ( KP != K ) {
               if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );
               cswap(KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA );

               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
               TEMP = A( K, K-1 );
               A[K, K-1] = A( KP, K-1 );
               A[KP, K-1] = TEMP;
            }

            K = K - 1;
            KP = -IPIV( K );
            if ( KP != K ) {
               if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );
               cswap(KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA );
               TEMP = A( K, K );
               A[K][K] = A( KP, KP );
               A[KP][KP] = TEMP;
            }
         }

         K = K - 1;
         GO TO 50;
         } // 60
      }

      return;
      }

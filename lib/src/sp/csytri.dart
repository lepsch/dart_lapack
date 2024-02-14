      void csytri(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Array<double> _WORK_, final Box<int> INFO,) {
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

      Complex            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
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
      // INTRINSIC ABS, MAX

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
         xerbla('CSYTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         for (INFO = N; INFO >= 1; INFO--) { // 10
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == ZERO ) return;
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) { // 20
            if( IPIV( INFO ) > 0 && A( INFO, INFO ) == ZERO ) return;
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

            A[K][K] = ONE / A( K, K );

            // Compute column K of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               csymv[UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K][K] = A( K, K ) - CDOTU( K-1, WORK, 1, A( 1, K ), 1 );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = A( K, K+1 );
            AK = A( K, K ) / T;
            AKP1 = A( K+1, K+1 ) / T;
            AKKP1 = A( K, K+1 ) / T;
            D = T*( AK*AKP1-ONE );
            A[K][K] = AKP1 / D;
            A[K+1][K+1] = AK / D;
            A[K][K+1] = -AKKP1 / D;

            // Compute columns K and K+1 of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               csymv[UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K][K] = A( K, K ) - CDOTU( K-1, WORK, 1, A( 1, K ), 1 )                A( K, K+1 ) = A( K, K+1 ) - CDOTU( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 );
               ccopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               csymv[UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K+1 ), 1 )                A( K+1, K+1] = A( K+1, K+1 ) - CDOTU( K-1, WORK, 1, A( 1, K+1 ), 1 );
            }
            KSTEP = 2;
         }

         KP = ( IPIV( K ) ).abs();
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the leading
            // submatrix A(1:k+1,1:k+1)

            cswap(KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
            cswap(K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA );
            TEMP = A( K, K );
            A[K][K] = A( KP, KP );
            A[KP][KP] = TEMP;
            if ( KSTEP == 2 ) {
               TEMP = A( K, K+1 );
               A[K][K+1] = A( KP, K+1 );
               A[KP][K+1] = TEMP;
            }
         }

         K = K + KSTEP;
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

            A[K][K] = ONE / A( K, K );

            // Compute column K of the inverse.

            if ( K < N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               csymv[UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K][K] = A( K, K ) - CDOTU( N-K, WORK, 1, A( K+1, K ), 1 );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = A( K, K-1 );
            AK = A( K-1, K-1 ) / T;
            AKP1 = A( K, K ) / T;
            AKKP1 = A( K, K-1 ) / T;
            D = T*( AK*AKP1-ONE );
            A[K-1][K-1] = AKP1 / D;
            A[K][K] = AK / D;
            A[K][K-1] = -AKKP1 / D;

            // Compute columns K-1 and K of the inverse.

            if ( K < N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               csymv[UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K][K] = A( K, K ) - CDOTU( N-K, WORK, 1, A( K+1, K ), 1 )                A( K, K-1 ) = A( K, K-1 ) - CDOTU( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               ccopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               csymv[UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1] = A( K-1, K-1 ) - CDOTU( N-K, WORK, 1, A( K+1, K-1 ), 1 );
            }
            KSTEP = 2;
         }

         KP = ( IPIV( K ) ).abs();
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the trailing
            // submatrix A(k-1:n,k-1:n)

            if (KP < N) cswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );
            cswap(KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA );
            TEMP = A( K, K );
            A[K][K] = A( KP, KP );
            A[KP][KP] = TEMP;
            if ( KSTEP == 2 ) {
               TEMP = A( K, K-1 );
               A[K][K-1] = A( KP, K-1 );
               A[KP][K-1] = TEMP;
            }
         }

         K = K - KSTEP;
         GO TO 50;
         } // 60
      }

      }

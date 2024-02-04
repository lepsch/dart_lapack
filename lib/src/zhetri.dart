      void zhetri(UPLO, N, A, LDA, IPIV, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex         A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      Complex         CONE, ZERO;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KP, KSTEP;
      double             AK, AKP1, D, T;
      Complex         AKKP1, TEMP;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- Complex         ZDOTC;
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZHEMV, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZHETRI', -INFO );
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

         // Compute inv(A) from the factorization A = U*D*U**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1;
         } // 30

         // If K > N, exit from loop.

         if (K > N) GO TO 50;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A[K, K] = ONE / (A( K, K )).toDouble();

            // Compute column K of the inverse.

            if ( K > 1 ) {
               zcopy(K-1, A( 1, K ), 1, WORK, 1 );
               zhemv[UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K] = A( K, K ) - DBLE( ZDOTC( K-1, WORK, 1, A( 1, K ), 1 ) );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ( A( K, K+1 ) ).abs();
            AK = (A( K, K )).toDouble() / T;
            AKP1 = (A( K+1, K+1 )).toDouble() / T;
            AKKP1 = A( K, K+1 ) / T;
            D = T*( AK*AKP1-ONE );
            A[K, K] = AKP1 / D;
            A[K+1, K+1] = AK / D;
            A[K, K+1] = -AKKP1 / D;

            // Compute columns K and K+1 of the inverse.

            if ( K > 1 ) {
               zcopy(K-1, A( 1, K ), 1, WORK, 1 );
               zhemv[UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K] = A( K, K ) - DBLE( ZDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )                A( K, K+1 ) = A( K, K+1 ) - ZDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 );
               zcopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               zhemv[UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO, A( 1, K+1 ), 1 )                A( K+1, K+1] = A( K+1, K+1 ) - DBLE( ZDOTC( K-1, WORK, 1, A( 1, K+1 ), 1 ) );
            }
            KSTEP = 2;
         }

         KP = ( IPIV( K ) ).abs();
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the leading
            // submatrix A(1:k+1,1:k+1)

            zswap(KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
            for (J = KP + 1; J <= K - 1; J++) { // 40
               TEMP = DCONJG( A( J, K ) );
               A[J, K] = DCONJG( A( KP, J ) );
               A[KP, J] = TEMP;
            } // 40
            A[KP, K] = DCONJG( A( KP, K ) );
            TEMP = A( K, K );
            A[K, K] = A( KP, KP );
            A[KP, KP] = TEMP;
            if ( KSTEP == 2 ) {
               TEMP = A( K, K+1 );
               A[K, K+1] = A( KP, K+1 );
               A[KP, K+1] = TEMP;
            }
         }

         K = K + KSTEP;
         GO TO 30;
         } // 50

      } else {

         // Compute inv(A) from the factorization A = L*D*L**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N;
         } // 60

         // If K < 1, exit from loop.

         if (K < 1) GO TO 80;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A[K, K] = ONE / (A( K, K )).toDouble();

            // Compute column K of the inverse.

            if ( K < N ) {
               zcopy(N-K, A( K+1, K ), 1, WORK, 1 );
               zhemv[UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K] = A( K, K ) - DBLE( ZDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ( A( K, K-1 ) ).abs();
            AK = (A( K-1, K-1 )).toDouble() / T;
            AKP1 = (A( K, K )).toDouble() / T;
            AKKP1 = A( K, K-1 ) / T;
            D = T*( AK*AKP1-ONE );
            A[K-1, K-1] = AKP1 / D;
            A[K, K] = AK / D;
            A[K, K-1] = -AKKP1 / D;

            // Compute columns K-1 and K of the inverse.

            if ( K < N ) {
               zcopy(N-K, A( K+1, K ), 1, WORK, 1 );
               zhemv[UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K] = A( K, K ) - DBLE( ZDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )                A( K, K-1 ) = A( K, K-1 ) - ZDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               zcopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               zhemv[UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1] = A( K-1, K-1 ) - DBLE( ZDOTC( N-K, WORK, 1, A( K+1, K-1 ), 1 ) );
            }
            KSTEP = 2;
         }

         KP = ( IPIV( K ) ).abs();
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the trailing
            // submatrix A(k-1:n,k-1:n)

            if (KP < N) zswap( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );
            for (J = K + 1; J <= KP - 1; J++) { // 70
               TEMP = DCONJG( A( J, K ) );
               A[J, K] = DCONJG( A( KP, J ) );
               A[KP, J] = TEMP;
            } // 70
            A[KP, K] = DCONJG( A( KP, K ) );
            TEMP = A( K, K );
            A[K, K] = A( KP, KP );
            A[KP, KP] = TEMP;
            if ( KSTEP == 2 ) {
               TEMP = A( K, K-1 );
               A[K, K-1] = A( KP, K-1 );
               A[KP, K-1] = TEMP;
            }
         }

         K = K - KSTEP;
         GO TO 60;
         } // 80
      }

      return;
      }

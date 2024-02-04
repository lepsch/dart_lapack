      void chptri(UPLO, N, AP, IPIV, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex            AP( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE;
      Complex            CONE, ZERO;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP;
      double               AK, AKP1, D, T;
      Complex            AKKP1, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHPMV, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('CHPTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         KP = N*( N+1 ) / 2;
         for (INFO = N; INFO >= 1; INFO--) { // 10
            if( IPIV( INFO ) > 0 && AP( KP ) == ZERO ) return;
            KP = KP - INFO;
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         KP = 1;
         for (INFO = 1; INFO <= N; INFO++) { // 20
            if( IPIV( INFO ) > 0 && AP( KP ) == ZERO ) return;
            KP = KP + N - INFO + 1;
         } // 20
      }
      INFO = 0;

      if ( UPPER ) {

         // Compute inv(A) from the factorization A = U*D*U**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1;
         KC = 1;
         } // 30

         // If K > N, exit from loop.

         if (K > N) GO TO 50;

         KCNEXT = KC + K;
         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            AP[KC+K-1] = ONE / REAL( AP( KC+K-1 ) );

            // Compute column K of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, AP( KC ), 1, WORK, 1 );
               chpmv[UPLO, K-1, -CONE, AP, WORK, 1, ZERO, AP( KC ), 1 )                AP( KC+K-1] = AP( KC+K-1 ) - double( CDOTC( K-1, WORK, 1, AP( KC ), 1 ) );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ( AP( KCNEXT+K-1 ) ).abs();
            AK = double( AP( KC+K-1 ) ) / T;
            AKP1 = double( AP( KCNEXT+K ) ) / T;
            AKKP1 = AP( KCNEXT+K-1 ) / T;
            D = T*( AK*AKP1-ONE );
            AP[KC+K-1] = AKP1 / D;
            AP[KCNEXT+K] = AK / D;
            AP[KCNEXT+K-1] = -AKKP1 / D;

            // Compute columns K and K+1 of the inverse.

            if ( K > 1 ) {
               ccopy(K-1, AP( KC ), 1, WORK, 1 );
               chpmv[UPLO, K-1, -CONE, AP, WORK, 1, ZERO, AP( KC ), 1 )                AP( KC+K-1] = AP( KC+K-1 ) - double( CDOTC( K-1, WORK, 1, AP( KC ), 1 ) )                AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) - CDOTC( K-1, AP( KC ), 1, AP( KCNEXT ), 1 );
               ccopy(K-1, AP( KCNEXT ), 1, WORK, 1 );
               chpmv[UPLO, K-1, -CONE, AP, WORK, 1, ZERO, AP( KCNEXT ), 1 )                AP( KCNEXT+K] = AP( KCNEXT+K ) - double( CDOTC( K-1, WORK, 1, AP( KCNEXT ), 1 ) );
            }
            KSTEP = 2;
            KCNEXT = KCNEXT + K + 1;
         }

         KP = ( IPIV( K ) ).abs();
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the leading
            // submatrix A(1:k+1,1:k+1)

            KPC = ( KP-1 )*KP / 2 + 1;
            cswap(KP-1, AP( KC ), 1, AP( KPC ), 1 );
            KX = KPC + KP - 1;
            for (J = KP + 1; J <= K - 1; J++) { // 40
               KX = KX + J - 1;
               TEMP = CONJG( AP( KC+J-1 ) );
               AP[KC+J-1] = CONJG( AP( KX ) );
               AP[KX] = TEMP;
            } // 40
            AP[KC+KP-1] = CONJG( AP( KC+KP-1 ) );
            TEMP = AP( KC+K-1 );
            AP[KC+K-1] = AP( KPC+KP-1 );
            AP[KPC+KP-1] = TEMP;
            if ( KSTEP == 2 ) {
               TEMP = AP( KC+K+K-1 );
               AP[KC+K+K-1] = AP( KC+K+KP-1 );
               AP[KC+K+KP-1] = TEMP;
            }
         }

         K = K + KSTEP;
         KC = KCNEXT;
         GO TO 30;
         } // 50

      } else {

         // Compute inv(A) from the factorization A = L*D*L**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         NPP = N*( N+1 ) / 2;
         K = N;
         KC = NPP;
         } // 60

         // If K < 1, exit from loop.

         if (K < 1) GO TO 80;

         KCNEXT = KC - ( N-K+2 );
         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            AP[KC] = ONE / REAL( AP( KC ) );

            // Compute column K of the inverse.

            if ( K < N ) {
               ccopy(N-K, AP( KC+1 ), 1, WORK, 1 );
               chpmv[UPLO, N-K, -CONE, AP( KC+N-K+1 ), WORK, 1, ZERO, AP( KC+1 ), 1 )                AP( KC] = AP( KC ) - double( CDOTC( N-K, WORK, 1, AP( KC+1 ), 1 ) );
            }
            KSTEP = 1;
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ( AP( KCNEXT+1 ) ).abs();
            AK = double( AP( KCNEXT ) ) / T;
            AKP1 = double( AP( KC ) ) / T;
            AKKP1 = AP( KCNEXT+1 ) / T;
            D = T*( AK*AKP1-ONE );
            AP[KCNEXT] = AKP1 / D;
            AP[KC] = AK / D;
            AP[KCNEXT+1] = -AKKP1 / D;

            // Compute columns K-1 and K of the inverse.

            if ( K < N ) {
               ccopy(N-K, AP( KC+1 ), 1, WORK, 1 );
               chpmv[UPLO, N-K, -CONE, AP( KC+( N-K+1 ) ), WORK, 1, ZERO, AP( KC+1 ), 1 )                AP( KC] = AP( KC ) - double( CDOTC( N-K, WORK, 1, AP( KC+1 ), 1 ) )                AP( KCNEXT+1 ) = AP( KCNEXT+1 ) - CDOTC( N-K, AP( KC+1 ), 1, AP( KCNEXT+2 ), 1 );
               ccopy(N-K, AP( KCNEXT+2 ), 1, WORK, 1 );
               chpmv[UPLO, N-K, -CONE, AP( KC+( N-K+1 ) ), WORK, 1, ZERO, AP( KCNEXT+2 ), 1 )                AP( KCNEXT] = AP( KCNEXT ) - double( CDOTC( N-K, WORK, 1, AP( KCNEXT+2 ), 1 ) );
            }
            KSTEP = 2;
            KCNEXT = KCNEXT - ( N-K+3 );
         }

         KP = ( IPIV( K ) ).abs();
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the trailing
            // submatrix A(k-1:n,k-1:n)

            KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1;
            if (KP < N) cswap( N-KP, AP( KC+KP-K+1 ), 1, AP( KPC+1 ), 1 );
            KX = KC + KP - K;
            for (J = K + 1; J <= KP - 1; J++) { // 70
               KX = KX + N - J + 1;
               TEMP = CONJG( AP( KC+J-K ) );
               AP[KC+J-K] = CONJG( AP( KX ) );
               AP[KX] = TEMP;
            } // 70
            AP[KC+KP-K] = CONJG( AP( KC+KP-K ) );
            TEMP = AP( KC );
            AP[KC] = AP( KPC );
            AP[KPC] = TEMP;
            if ( KSTEP == 2 ) {
               TEMP = AP( KC-N+K-1 );
               AP[KC-N+K-1] = AP( KC-N+KP-1 );
               AP[KC-N+KP-1] = TEMP;
            }
         }

         K = K - KSTEP;
         KC = KCNEXT;
         GO TO 60;
         } // 80
      }

      return;
      }
      SUBROUTINE ZHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      COMPLEX*16         CONE, CZERO
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KP, KSTEP;
      double             AK, AKP1, D, T;
      COMPLEX*16         AKKP1, TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZHEMV, ZSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, DBLE
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('ZHETRI_ROOK', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ) > 0 && A( INFO, INFO ) == CZERO ) RETURN
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) { // 20
            IF( IPIV( INFO ) > 0 && A( INFO, INFO ) == CZERO ) RETURN
         } // 20
      }
      INFO = 0

      if ( UPPER ) {

         // Compute inv(A) from the factorization A = U*D*U**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1
         } // 30

         // If K > N, exit from loop.

         if (K > N) GO TO 70;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / DBLE( A( K, K ) )

            // Compute column K of the inverse.

            if ( K > 1 ) {
               zcopy(K-1, A( 1, K ), 1, WORK, 1 );
               zhemv(UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - DBLE( ZDOTC( K-1, WORK, 1, A( 1, K ), 1 ) );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K+1 ) )
            AK = DBLE( A( K, K ) ) / T
            AKP1 = DBLE( A( K+1, K+1 ) ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D

            // Compute columns K and K+1 of the inverse.

            if ( K > 1 ) {
               zcopy(K-1, A( 1, K ), 1, WORK, 1 );
               zhemv(UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K ), 1 );
               A( K, K ) = A( K, K ) - DBLE( ZDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )                A( K, K+1 ) = A( K, K+1 ) - ZDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               zcopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               zhemv(UPLO, K-1, -CONE, A, LDA, WORK, 1, CZERO, A( 1, K+1 ), 1 )                A( K+1, K+1 ) = A( K+1, K+1 ) - DBLE( ZDOTC( K-1, WORK, 1, A( 1, K+1 ), 1 ) );
            }
            KSTEP = 2
         }

         if ( KSTEP == 1 ) {

            // Interchange rows and columns K and IPIV(K) in the leading
            // submatrix A(1:k,1:k)

            KP = IPIV( K )
            if ( KP != K ) {

               if (KP > 1) CALL ZSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );

               for (J = KP + 1; J <= K - 1; J++) { // 40
                  TEMP = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
               } // 40

               A( KP, K ) = DCONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            }
         } else {

            // Interchange rows and columns K and K+1 with -IPIV(K) and
            // -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n)

            // (1) Interchange rows and columns K and -IPIV(K)

            KP = -IPIV( K )
            if ( KP != K ) {

               if (KP > 1) CALL ZSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );

               for (J = KP + 1; J <= K - 1; J++) { // 50
                  TEMP = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
               } // 50

               A( KP, K ) = DCONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP

               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            }

            // (2) Interchange rows and columns K+1 and -IPIV(K+1)

            K = K + 1
            KP = -IPIV( K )
            if ( KP != K ) {

               if (KP > 1) CALL ZSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 );

               for (J = KP + 1; J <= K - 1; J++) { // 60
                  TEMP = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
               } // 60

               A( KP, K ) = DCONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            }
         }

         K = K + 1
         GO TO 30
         } // 70

      } else {

         // Compute inv(A) from the factorization A = L*D*L**H.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N
         } // 80

         // If K < 1, exit from loop.

         if (K < 1) GO TO 120;

         if ( IPIV( K ) > 0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / DBLE( A( K, K ) )

            // Compute column K of the inverse.

            if ( K < N ) {
               zcopy(N-K, A( K+1, K ), 1, WORK, 1 );
               zhemv(UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - DBLE( ZDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K-1 ) )
            AK = DBLE( A( K-1, K-1 ) ) / T
            AKP1 = DBLE( A( K, K ) ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D

            // Compute columns K-1 and K of the inverse.

            if ( K < N ) {
               zcopy(N-K, A( K+1, K ), 1, WORK, 1 );
               zhemv(UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - DBLE( ZDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )                A( K, K-1 ) = A( K, K-1 ) - ZDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               zcopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               zhemv(UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, CZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1 ) = A( K-1, K-1 ) - DBLE( ZDOTC( N-K, WORK, 1, A( K+1, K-1 ), 1 ) );
            }
            KSTEP = 2
         }

         if ( KSTEP == 1 ) {

            // Interchange rows and columns K and IPIV(K) in the trailing
            // submatrix A(k:n,k:n)

            KP = IPIV( K )
            if ( KP != K ) {

               if (KP < N) CALL ZSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );

               for (J = K + 1; J <= KP - 1; J++) { // 90
                  TEMP = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
               } // 90

               A( KP, K ) = DCONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            }
         } else {

            // Interchange rows and columns K and K-1 with -IPIV(K) and
            // -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n)

            // (1) Interchange rows and columns K and -IPIV(K)

            KP = -IPIV( K )
            if ( KP != K ) {

               if (KP < N) CALL ZSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );

               for (J = K + 1; J <= KP - 1; J++) { // 100
                  TEMP = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
              } // 100

               A( KP, K ) = DCONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP

               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            }

            // (2) Interchange rows and columns K-1 and -IPIV(K-1)

            K = K - 1
            KP = -IPIV( K )
            if ( KP != K ) {

               if (KP < N) CALL ZSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );

               for (J = K + 1; J <= KP - 1; J++) { // 110
                  TEMP = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = TEMP
              } // 110

               A( KP, K ) = DCONJG( A( KP, K ) )

               TEMP = A( K, K )
               A( K, K ) = A( KP, KP )
               A( KP, KP ) = TEMP
            }
         }

         K = K - 1
         GO TO 80
         } // 120
      }

      RETURN

      // End of ZHETRI_ROOK

      }

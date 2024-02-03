      SUBROUTINE SSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                K, KP, KSTEP;
      REAL               AK, AKKP1, AKP1, D, T, TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT
      // EXTERNAL LSAME, SDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SSWAP, SSYMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('SSYTRI', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 && A( INFO, INFO ) == ZERO ) RETURN
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) { // 20
            IF( IPIV( INFO ).GT.0 && A( INFO, INFO ) == ZERO ) RETURN
         } // 20
      }
      INFO = 0

      if ( UPPER ) {

         // Compute inv(A) from the factorization A = U*D*U**T.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1
         } // 30

         // If K > N, exit from loop.

         if (K.GT.N) GO TO 40;

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / A( K, K )

            // Compute column K of the inverse.

            if ( K.GT.1 ) {
               scopy(K-1, A( 1, K ), 1, WORK, 1 );
               ssymv(UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( K-1, WORK, 1, A( 1, K ), 1 );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K+1 ) )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D

            // Compute columns K and K+1 of the inverse.

            if ( K.GT.1 ) {
               scopy(K-1, A( 1, K ), 1, WORK, 1 );
               ssymv(UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( K-1, WORK, 1, A( 1, K ), 1 )                A( K, K+1 ) = A( K, K+1 ) - SDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 );
               scopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               ssymv(UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO, A( 1, K+1 ), 1 )                A( K+1, K+1 ) = A( K+1, K+1 ) - SDOT( K-1, WORK, 1, A( 1, K+1 ), 1 );
            }
            KSTEP = 2
         }

         KP = ABS( IPIV( K ) )
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the leading
            // submatrix A(1:k+1,1:k+1)

            sswap(KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
            sswap(K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA );
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            if ( KSTEP == 2 ) {
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            }
         }

         K = K + KSTEP
         GO TO 30
         } // 40

      } else {

         // Compute inv(A) from the factorization A = L*D*L**T.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N
         } // 50

         // If K < 1, exit from loop.

         if (K.LT.1) GO TO 60;

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / A( K, K )

            // Compute column K of the inverse.

            if ( K.LT.N ) {
               scopy(N-K, A( K+1, K ), 1, WORK, 1 );
               ssymv(UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( N-K, WORK, 1, A( K+1, K ), 1 );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D

            // Compute columns K-1 and K of the inverse.

            if ( K.LT.N ) {
               scopy(N-K, A( K+1, K ), 1, WORK, 1 );
               ssymv(UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - SDOT( N-K, WORK, 1, A( K+1, K ), 1 )                A( K, K-1 ) = A( K, K-1 ) - SDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               scopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               ssymv(UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1 ) = A( K-1, K-1 ) - SDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 );
            }
            KSTEP = 2
         }

         KP = ABS( IPIV( K ) )
         if ( KP != K ) {

            // Interchange rows and columns K and KP in the trailing
            // submatrix A(k-1:n,k-1:n)

            if (KP.LT.N) CALL SSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 );
            sswap(KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA );
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            if ( KSTEP == 2 ) {
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            }
         }

         K = K - KSTEP
         GO TO 50
         } // 60
      }

      RETURN

      // End of SSYTRI

      }

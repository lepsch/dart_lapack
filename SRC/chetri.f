      SUBROUTINE CHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      COMPLEX            CONE, ZERO
      const              ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KP, KSTEP;
      REAL               AK, AKP1, D, T
      COMPLEX            AKKP1, TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHEMV, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('CHETRI', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (INFO = 1; INFO <= N; INFO++) { // 20
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO ) RETURN
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

         IF( K.GT.N ) GO TO 50

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / REAL( A( K, K ) )

            // Compute column K of the inverse.

            if ( K.GT.1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               chemv(UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K+1 ) )
            AK = REAL( A( K, K ) ) / T
            AKP1 = REAL( A( K+1, K+1 ) ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D

            // Compute columns K and K+1 of the inverse.

            if ( K.GT.1 ) {
               ccopy(K-1, A( 1, K ), 1, WORK, 1 );
               chemv(UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO, A( 1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K ), 1 ) )                A( K, K+1 ) = A( K, K+1 ) - CDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 );
               ccopy(K-1, A( 1, K+1 ), 1, WORK, 1 );
               chemv(UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO, A( 1, K+1 ), 1 )                A( K+1, K+1 ) = A( K+1, K+1 ) - REAL( CDOTC( K-1, WORK, 1, A( 1, K+1 ), 1 ) );
            }
            KSTEP = 2
         }

         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) {

            // Interchange rows and columns K and KP in the leading
            // submatrix A(1:k+1,1:k+1)

            cswap(KP-1, A( 1, K ), 1, A( 1, KP ), 1 );
            for (J = KP + 1; J <= K - 1; J++) { // 40
               TEMP = CONJG( A( J, K ) )
               A( J, K ) = CONJG( A( KP, J ) )
               A( KP, J ) = TEMP
            } // 40
            A( KP, K ) = CONJG( A( KP, K ) )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            if ( KSTEP.EQ.2 ) {
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            }
         }

         K = K + KSTEP
         GO TO 30
         } // 50

      } else {

         // Compute inv(A) from the factorization A = L*D*L**H.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N
         } // 60

         // If K < 1, exit from loop.

         IF( K.LT.1 ) GO TO 80

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Invert the diagonal block.

            A( K, K ) = ONE / REAL( A( K, K ) )

            // Compute column K of the inverse.

            if ( K.LT.N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               chemv(UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) );
            }
            KSTEP = 1
         } else {

            // 2 x 2 diagonal block

            // Invert the diagonal block.

            T = ABS( A( K, K-1 ) )
            AK = REAL( A( K-1, K-1 ) ) / T
            AKP1 = REAL( A( K, K ) ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D

            // Compute columns K-1 and K of the inverse.

            if ( K.LT.N ) {
               ccopy(N-K, A( K+1, K ), 1, WORK, 1 );
               chemv(UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K ), 1 )                A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K ), 1 ) )                A( K, K-1 ) = A( K, K-1 ) - CDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 );
               ccopy(N-K, A( K+1, K-1 ), 1, WORK, 1 );
               chemv(UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK, 1, ZERO, A( K+1, K-1 ), 1 )                A( K-1, K-1 ) = A( K-1, K-1 ) - REAL( CDOTC( N-K, WORK, 1, A( K+1, K-1 ), 1 ) );
            }
            KSTEP = 2
         }

         KP = ABS( IPIV( K ) )
         if ( KP.NE.K ) {

            // Interchange rows and columns K and KP in the trailing
            // submatrix A(k-1:n,k-1:n)

            IF( KP.LT.N ) CALL CSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            for (J = K + 1; J <= KP - 1; J++) { // 70
               TEMP = CONJG( A( J, K ) )
               A( J, K ) = CONJG( A( KP, J ) )
               A( KP, J ) = TEMP
            } // 70
            A( KP, K ) = CONJG( A( KP, K ) )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            if ( KSTEP.EQ.2 ) {
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            }
         }

         K = K - KSTEP
         GO TO 60
         } // 80
      }

      RETURN

      // End of CHETRI

      }

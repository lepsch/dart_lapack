      SUBROUTINE CHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, K, KP;
      REAL               S
      COMPLEX            AK, AKM1, AKM1K, BK, BKM1, DENOM
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERU, CLACGV, CSSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, REAL
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CHETRS', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      if ( UPPER ) {

         // Solve A*X = B, where A = U*D*U**H.

         // First solve U*D*X = B, overwriting B with X.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N
   10    CONTINUE

         // If K < 1, exit from loop.

         IF( K.LT.1 ) GO TO 30

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Interchange rows K and IPIV(K).

            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )

            // Multiply by inv(U(K)), where U(K) is the transformation
            // stored in column K of A.

            CALL CGERU( K-1, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB )

            // Multiply by the inverse of the diagonal block.

            S = REAL( ONE ) / REAL( A( K, K ) )
            CALL CSSCAL( NRHS, S, B( K, 1 ), LDB )
            K = K - 1
         } else {

            // 2 x 2 diagonal block

            // Interchange rows K-1 and -IPIV(K).

            KP = -IPIV( K )
            IF( KP.NE.K-1 ) CALL CSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )

            // Multiply by inv(U(K)), where U(K) is the transformation
            // stored in columns K-1 and K of A.

            CALL CGERU( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB, B( 1, 1 ), LDB )             CALL CGERU( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ), LDB, B( 1, 1 ), LDB )

            // Multiply by the inverse of the diagonal block.

            AKM1K = A( K-1, K )
            AKM1 = A( K-1, K-1 ) / AKM1K
            AK = A( K, K ) / CONJG( AKM1K )
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / CONJG( AKM1K )
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            K = K - 2
         }

         GO TO 10
   30    CONTINUE

         // Next solve U**H *X = B, overwriting B with X.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1
   40    CONTINUE

         // If K > N, exit from loop.

         IF( K.GT.N ) GO TO 50

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Multiply by inv(U**H(K)), where U(K) is the transformation
            // stored in column K of A.

            if ( K.GT.1 ) {
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
            }

            // Interchange rows K and IPIV(K).

            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 1
         } else {

            // 2 x 2 diagonal block

            // Multiply by inv(U**H(K+1)), where U(K+1) is the transformation
            // stored in columns K and K+1 of A.

            if ( K.GT.1 ) {
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ), 1, ONE, B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )

               CALL CLACGV( NRHS, B( K+1, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
               CALL CLACGV( NRHS, B( K+1, 1 ), LDB )
            }

            // Interchange rows K and -IPIV(K).

            KP = -IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 2
         }

         GO TO 40
   50    CONTINUE

      } else {

         // Solve A*X = B, where A = L*D*L**H.

         // First solve L*D*X = B, overwriting B with X.

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = 1
   60    CONTINUE

         // If K > N, exit from loop.

         IF( K.GT.N ) GO TO 80

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Interchange rows K and IPIV(K).

            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )

            // Multiply by inv(L(K)), where L(K) is the transformation
            // stored in column K of A.

            IF( K.LT.N ) CALL CGERU( N-K, NRHS, -ONE, A( K+1, K ), 1, B( K, 1 ), LDB, B( K+1, 1 ), LDB )

            // Multiply by the inverse of the diagonal block.

            S = REAL( ONE ) / REAL( A( K, K ) )
            CALL CSSCAL( NRHS, S, B( K, 1 ), LDB )
            K = K + 1
         } else {

            // 2 x 2 diagonal block

            // Interchange rows K+1 and -IPIV(K).

            KP = -IPIV( K )
            IF( KP.NE.K+1 ) CALL CSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )

            // Multiply by inv(L(K)), where L(K) is the transformation
            // stored in columns K and K+1 of A.

            if ( K.LT.N-1 ) {
               CALL CGERU( N-K-1, NRHS, -ONE, A( K+2, K ), 1, B( K, 1 ), LDB, B( K+2, 1 ), LDB )                CALL CGERU( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1, B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            }

            // Multiply by the inverse of the diagonal block.

            AKM1K = A( K+1, K )
            AKM1 = A( K, K ) / CONJG( AKM1K )
            AK = A( K+1, K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / CONJG( AKM1K )
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            K = K + 2
         }

         GO TO 60
   80    CONTINUE

         // Next solve L**H *X = B, overwriting B with X.

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2, depending on the size of the diagonal blocks.

         K = N
   90    CONTINUE

         // If K < 1, exit from loop.

         IF( K.LT.1 ) GO TO 100

         if ( IPIV( K ).GT.0 ) {

            // 1 x 1 diagonal block

            // Multiply by inv(L**H(K)), where L(K) is the transformation
            // stored in column K of A.

            if ( K.LT.N ) {
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
            }

            // Interchange rows K and IPIV(K).

            KP = IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         } else {

            // 2 x 2 diagonal block

            // Multiply by inv(L**H(K-1)), where L(K-1) is the transformation
            // stored in columns K-1 and K of A.

            if ( K.LT.N ) {
               CALL CLACGV( NRHS, B( K, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
               CALL CLACGV( NRHS, B( K, 1 ), LDB )

               CALL CLACGV( NRHS, B( K-1, 1 ), LDB )
               CALL CGEMV( 'Conjugate transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, A( K+1, K-1 ), 1, ONE, B( K-1, 1 ), LDB )
               CALL CLACGV( NRHS, B( K-1, 1 ), LDB )
            }

            // Interchange rows K and -IPIV(K).

            KP = -IPIV( K )
            IF( KP.NE.K ) CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 2
         }

         GO TO 90
  100    CONTINUE
      }

      RETURN

      // End of CHETRS

      }

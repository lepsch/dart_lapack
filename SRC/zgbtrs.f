      SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AB( LDAB, * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LNOTI, NOTRAN;
      int                I, J, KD, L, LM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZGERU, ZLACGV, ZSWAP, ZTBSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NOTRAN = LSAME( TRANS, 'N' );
      if ( !NOTRAN && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 ) {
         INFO = -3;
      } else if ( KU < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDAB < ( 2*KL+KU+1 ) ) {
         INFO = -7;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('ZGBTRS', -INFO );
         RETURN;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      KD = KU + KL + 1;
      LNOTI = KL > 0;

      if ( NOTRAN ) {

         // Solve  A*X = B.

         // Solve L*X = B, overwriting B with X.

         // L is represented as a product of permutations and unit lower
         // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
         // where each transformation L(i) is a rank-one modification of
         // the identity matrix.

         if ( LNOTI ) {
            for (J = 1; J <= N - 1; J++) { // 10
               LM = MIN( KL, N-J );
               L = IPIV( J );
               if (L != J) CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB );
               zgeru(LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), LDB, B( J+1, 1 ), LDB );
            } // 10
         }

         for (I = 1; I <= NRHS; I++) { // 20

            // Solve U*X = B, overwriting B with X.

            ztbsv('Upper', 'No transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
         } // 20

      } else if ( LSAME( TRANS, 'T' ) ) {

         // Solve A**T * X = B.

         for (I = 1; I <= NRHS; I++) { // 30

            // Solve U**T * X = B, overwriting B with X.

            ztbsv('Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
         } // 30

         // Solve L**T * X = B, overwriting B with X.

         if ( LNOTI ) {
            DO 40 J = N - 1, 1, -1;
               LM = MIN( KL, N-J );
               zgemv('Transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB );
               L = IPIV( J );
               if (L != J) CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB );
            } // 40
         }

      } else {

         // Solve A**H * X = B.

         for (I = 1; I <= NRHS; I++) { // 50

            // Solve U**H * X = B, overwriting B with X.

            ztbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
         } // 50

         // Solve L**H * X = B, overwriting B with X.

         if ( LNOTI ) {
            DO 60 J = N - 1, 1, -1;
               LM = MIN( KL, N-J );
               zlacgv(NRHS, B( J, 1 ), LDB );
               zgemv('Conjugate transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB );
               zlacgv(NRHS, B( J, 1 ), LDB );
               L = IPIV( J );
               if (L != J) CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB );
            } // 60
         }
      }
      RETURN;

      // End of ZGBTRS

      }

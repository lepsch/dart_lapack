      SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AB( LDAB, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
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

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 ) {
         INFO = -3
      } else if ( KU.LT.0 ) {
         INFO = -4
      } else if ( NRHS.LT.0 ) {
         INFO = -5
      } else if ( LDAB.LT.( 2*KL+KU+1 ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGBTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      KD = KU + KL + 1
      LNOTI = KL.GT.0

      if ( NOTRAN ) {

         // Solve  A*X = B.

         // Solve L*X = B, overwriting B with X.

         // L is represented as a product of permutations and unit lower
         // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
         // where each transformation L(i) is a rank-one modification of
         // the identity matrix.

         if ( LNOTI ) {
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )                CALL ZGERU( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         }

         DO 20 I = 1, NRHS

            // Solve U*X = B, overwriting B with X.

            ztbsv('Upper', 'No transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
   20    CONTINUE

      } else if ( LSAME( TRANS, 'T' ) ) {

         // Solve A**T * X = B.

         DO 30 I = 1, NRHS

            // Solve U**T * X = B, overwriting B with X.

            ztbsv('Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
   30    CONTINUE

         // Solve L**T * X = B, overwriting B with X.

         if ( LNOTI ) {
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               zgemv('Transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB );
               L = IPIV( J )
               IF( L.NE.J ) CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         }

      } else {

         // Solve A**H * X = B.

         DO 50 I = 1, NRHS

            // Solve U**H * X = B, overwriting B with X.

            ztbsv('Upper', 'Conjugate transpose', 'Non-unit', N, KL+KU, AB, LDAB, B( 1, I ), 1 );
   50    CONTINUE

         // Solve L**H * X = B, overwriting B with X.

         if ( LNOTI ) {
            DO 60 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               zlacgv(NRHS, B( J, 1 ), LDB );
               zgemv('Conjugate transpose', LM, NRHS, -ONE, B( J+1, 1 ), LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB );
               zlacgv(NRHS, B( J, 1 ), LDB );
               L = IPIV( J )
               IF( L.NE.J ) CALL ZSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   60       CONTINUE
         }
      }
      RETURN

      // End of ZGBTRS

      }

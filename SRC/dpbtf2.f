      SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, KLD, KN;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('DPBTF2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      KLD = MAX( 1, LDAB-1 )

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**T*U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) GO TO 30
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ

            // Compute elements J+1:J+KN of row J and update the
            // trailing submatrix within the band.

            KN = MIN( KD, N-J )
            if ( KN.GT.0 ) {
               dscal(KN, ONE / AJJ, AB( KD, J+1 ), KLD );
               dsyr('Upper', KN, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD );
            }
   10    CONTINUE
      } else {

         // Compute the Cholesky factorization A = L*L**T.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) GO TO 30
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ

            // Compute elements J+1:J+KN of column J and update the
            // trailing submatrix within the band.

            KN = MIN( KD, N-J )
            if ( KN.GT.0 ) {
               dscal(KN, ONE / AJJ, AB( 2, J ), 1 );
               dsyr('Lower', KN, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD );
            }
   20    CONTINUE
      }
      RETURN

   30 CONTINUE
      INFO = J
      RETURN

      // End of DPBTF2

      }

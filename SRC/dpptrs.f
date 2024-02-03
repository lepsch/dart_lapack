      SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             AP( * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DPPTRS', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      if ( UPPER ) {

         // Solve A*X = B where A = U**T * U.

         DO 10 I = 1, NRHS

            // Solve U**T *X = B, overwriting B with X.

            CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 )

            // Solve U*X = B, overwriting B with X.

            CALL DTPSV( 'Upper', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
   10    CONTINUE
      } else {

         // Solve A*X = B where A = L * L**T.

         DO 20 I = 1, NRHS

            // Solve L*Y = B, overwriting B with X.

            CALL DTPSV( 'Lower', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 )

            // Solve L**T *X = Y, overwriting B with X.

            CALL DTPSV( 'Lower', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
   20    CONTINUE
      }

      RETURN

      // End of DPPTRS

      }

      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASWP, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( !NOTRAN && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('DGETRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      if ( NOTRAN ) {

         // Solve A * X = B.

         // Apply row interchanges to the right hand sides.

         dlaswp(NRHS, B, LDB, 1, N, IPIV, 1 );

         // Solve L*X = B, overwriting B with X.

         dtrsm('Left', 'Lower', 'No transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve U*X = B, overwriting B with X.

         dtrsm('Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );
      } else {

         // Solve A**T * X = B.

         // Solve U**T *X = B, overwriting B with X.

         dtrsm('Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Solve L**T *X = B, overwriting B with X.

         dtrsm('Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, A, LDA, B, LDB );

         // Apply row interchanges to the solution vectors.

         dlaswp(NRHS, B, LDB, 1, N, IPIV, -1 );
      }

      RETURN

      // End of DGETRS

      }

      SUBROUTINE ZTPTRI( UPLO, DIAG, N, AP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      const              ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JC, JCLAST, JJ;
      COMPLEX*16         AJJ
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSCAL, ZTPMV
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      }
      if ( INFO.NE.0 ) {
         xerbla('ZTPTRI', -INFO );
         RETURN
      }

      // Check for singularity if non-unit.

      if ( NOUNIT ) {
         if ( UPPER ) {
            JJ = 0
            for (INFO = 1; INFO <= N; INFO++) { // 10
               JJ = JJ + INFO
               IF( AP( JJ ).EQ.ZERO ) RETURN
   10       CONTINUE
         } else {
            JJ = 1
            for (INFO = 1; INFO <= N; INFO++) { // 20
               IF( AP( JJ ).EQ.ZERO ) RETURN
               JJ = JJ + N - INFO + 1
   20       CONTINUE
         }
         INFO = 0
      }

      if ( UPPER ) {

         // Compute inverse of upper triangular matrix.

         JC = 1
         for (J = 1; J <= N; J++) { // 30
            if ( NOUNIT ) {
               AP( JC+J-1 ) = ONE / AP( JC+J-1 )
               AJJ = -AP( JC+J-1 )
            } else {
               AJJ = -ONE
            }

            // Compute elements 1:j-1 of j-th column.

            ztpmv('Upper', 'No transpose', DIAG, J-1, AP, AP( JC ), 1 );
            zscal(J-1, AJJ, AP( JC ), 1 );
            JC = JC + J
   30    CONTINUE

      } else {

         // Compute inverse of lower triangular matrix.

         JC = N*( N+1 ) / 2
         DO 40 J = N, 1, -1
            if ( NOUNIT ) {
               AP( JC ) = ONE / AP( JC )
               AJJ = -AP( JC )
            } else {
               AJJ = -ONE
            }
            if ( J.LT.N ) {

               // Compute elements j+1:n of j-th column.

               ztpmv('Lower', 'No transpose', DIAG, N-J, AP( JCLAST ), AP( JC+1 ), 1 );
               zscal(N-J, AJJ, AP( JC+1 ), 1 );
            }
            JCLAST = JC
            JC = JC - N + J - 2
   40    CONTINUE
      }

      RETURN

      // End of ZTPTRI

      }

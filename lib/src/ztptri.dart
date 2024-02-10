      void ztptri(UPLO, DIAG, N, AP, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                INFO, N;
      Complex         AP( * );
      // ..

      Complex         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      bool               NOUNIT, UPPER;
      int                J, JC, JCLAST, JJ;
      Complex         AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSCAL, ZTPMV

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOUNIT = lsame( DIAG, 'N' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('ZTPTRI', -INFO );
         return;
      }

      // Check for singularity if non-unit.

      if ( NOUNIT ) {
         if ( UPPER ) {
            JJ = 0;
            for (INFO = 1; INFO <= N; INFO++) { // 10
               JJ = JJ + INFO;
               if( AP( JJ ) == ZERO ) return;
            } // 10
         } else {
            JJ = 1;
            for (INFO = 1; INFO <= N; INFO++) { // 20
               if( AP( JJ ) == ZERO ) return;
               JJ = JJ + N - INFO + 1;
            } // 20
         }
         INFO = 0;
      }

      if ( UPPER ) {

         // Compute inverse of upper triangular matrix.

         JC = 1;
         for (J = 1; J <= N; J++) { // 30
            if ( NOUNIT ) {
               AP[JC+J-1] = ONE / AP( JC+J-1 );
               AJJ = -AP( JC+J-1 );
            } else {
               AJJ = -ONE;
            }

            // Compute elements 1:j-1 of j-th column.

            ztpmv('Upper', 'No transpose', DIAG, J-1, AP, AP( JC ), 1 );
            zscal(J-1, AJJ, AP( JC ), 1 );
            JC = JC + J;
         } // 30

      } else {

         // Compute inverse of lower triangular matrix.

         JC = N*( N+1 ) / 2;
         for (J = N; J >= 1; J--) { // 40
            if ( NOUNIT ) {
               AP[JC] = ONE / AP( JC );
               AJJ = -AP( JC );
            } else {
               AJJ = -ONE;
            }
            if ( J < N ) {

               // Compute elements j+1:n of j-th column.

               ztpmv('Lower', 'No transpose', DIAG, N-J, AP( JCLAST ), AP( JC+1 ), 1 );
               zscal(N-J, AJJ, AP( JC+1 ), 1 );
            }
            JCLAST = JC;
            JC = JC - N + J - 2;
         } // 40
      }

      }

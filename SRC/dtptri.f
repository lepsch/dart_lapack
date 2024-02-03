      SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JC, JCLAST, JJ;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DTPMV, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      NOUNIT = LSAME( DIAG, 'N' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('DTPTRI', -INFO );
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
               AP( JC+J-1 ) = ONE / AP( JC+J-1 );
               AJJ = -AP( JC+J-1 );
            } else {
               AJJ = -ONE;
            }

            // Compute elements 1:j-1 of j-th column.

            dtpmv('Upper', 'No transpose', DIAG, J-1, AP, AP( JC ), 1 );
            dscal(J-1, AJJ, AP( JC ), 1 );
            JC = JC + J;
         } // 30

      } else {

         // Compute inverse of lower triangular matrix.

         JC = N*( N+1 ) / 2;
         DO 40 J = N, 1, -1;
            if ( NOUNIT ) {
               AP( JC ) = ONE / AP( JC );
               AJJ = -AP( JC );
            } else {
               AJJ = -ONE;
            }
            if ( J < N ) {

               // Compute elements j+1:n of j-th column.

               dtpmv('Lower', 'No transpose', DIAG, N-J, AP( JCLAST ), AP( JC+1 ), 1 );
               dscal(N-J, AJJ, AP( JC+1 ), 1 );
            }
            JCLAST = JC;
            JC = JC - N + J - 2;
         } // 40
      }

      return;

      // End of DTPTRI

      }

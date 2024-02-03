      void zpptrf(UPLO, N, AP, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      Complex         AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      Complex         ZDOTC;
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZHPR, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('ZPPTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**H * U.

         JJ = 0;
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1;
            JJ = JJ + J;

            // Compute elements 1:J-1 of column J.

            if (J > 1) ztpsv( 'Upper', 'Conjugate transpose', 'Non-unit', J-1, AP, AP( JC ), 1 );

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = DBLE( AP( JJ ) ) - DBLE( ZDOTC( J-1, AP( JC ), 1, AP( JC ), 1 ) );
            if ( AJJ <= ZERO ) {
               AP( JJ ) = AJJ;
               GO TO 30;
            }
            AP( JJ ) = sqrt( AJJ );
         } // 10
      } else {

         // Compute the Cholesky factorization A = L * L**H.

         JJ = 1;
         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = DBLE( AP( JJ ) );
            if ( AJJ <= ZERO ) {
               AP( JJ ) = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            AP( JJ ) = AJJ;

            // Compute elements J+1:N of column J and update the trailing
            // submatrix.

            if ( J < N ) {
               zdscal(N-J, ONE / AJJ, AP( JJ+1 ), 1 );
               zhpr('Lower', N-J, -ONE, AP( JJ+1 ), 1, AP( JJ+N-J+1 ) );
               JJ = JJ + N - J + 1;
            }
         } // 20
      }
      GO TO 40;

      } // 30
      INFO = J;

      } // 40
      return;
      }

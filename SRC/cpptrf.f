      SUBROUTINE CPPTRF( UPLO, N, AP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ;
      REAL               AJJ
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPR, CSSCAL, CTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('CPPTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**H * U.

         JJ = 0
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1
            JJ = JJ + J

            // Compute elements 1:J-1 of column J.

            if (J.GT.1) CALL CTPSV( 'Upper', 'Conjugate transpose', 'Non-unit', J-1, AP, AP( JC ), 1 );

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = REAL( REAL( AP( JJ ) ) - CDOTC( J-1, AP( JC ), 1, AP( JC ), 1 ) )
            if ( AJJ.LE.ZERO ) {
               AP( JJ ) = AJJ
               GO TO 30
            }
            AP( JJ ) = SQRT( AJJ )
         } // 10
      } else {

         // Compute the Cholesky factorization A = L * L**H.

         JJ = 1
         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = REAL( AP( JJ ) )
            if ( AJJ.LE.ZERO ) {
               AP( JJ ) = AJJ
               GO TO 30
            }
            AJJ = SQRT( AJJ )
            AP( JJ ) = AJJ

            // Compute elements J+1:N of column J and update the trailing
            // submatrix.

            if ( J.LT.N ) {
               csscal(N-J, ONE / AJJ, AP( JJ+1 ), 1 );
               chpr('Lower', N-J, -ONE, AP( JJ+1 ), 1, AP( JJ+N-J+1 ) );
               JJ = JJ + N - J + 1
            }
         } // 20
      }
      GO TO 40

      } // 30
      INFO = J

      } // 40
      RETURN

      // End of CPPTRF

      }

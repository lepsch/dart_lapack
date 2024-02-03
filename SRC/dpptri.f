      SUBROUTINE DPPTRI( UPLO, N, AP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             AP( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ, JJN;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT;
      // EXTERNAL LSAME, DDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSPR, DTPMV, DTPTRI, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      }
      if ( INFO != 0 ) {
         xerbla('DPPTRI', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Invert the triangular Cholesky factor U or L.

      dtptri(UPLO, 'Non-unit', N, AP, INFO );
      if (INFO > 0) RETURN;

      if ( UPPER ) {

         // Compute the product inv(U) * inv(U)**T.

         JJ = 0
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1
            JJ = JJ + J
            if (J > 1) CALL DSPR( 'Upper', J-1, ONE, AP( JC ), 1, AP );
            AJJ = AP( JJ )
            dscal(J, AJJ, AP( JC ), 1 );
         } // 10

      } else {

         // Compute the product inv(L)**T * inv(L).

         JJ = 1
         for (J = 1; J <= N; J++) { // 20
            JJN = JJ + N - J + 1
            AP( JJ ) = DDOT( N-J+1, AP( JJ ), 1, AP( JJ ), 1 )
            if (J < N) CALL DTPMV( 'Lower', 'Transpose', 'Non-unit', N-J, AP( JJN ), AP( JJ+1 ), 1 );
            JJ = JJN
         } // 20
      }

      RETURN

      // End of DPPTRI

      }

      SUBROUTINE SPPTRI( UPLO, N, AP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ, JJN;
      REAL               AJJ
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT
      // EXTERNAL LSAME, SDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSPR, STPMV, STPTRI, XERBLA
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
         xerbla('SPPTRI', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Invert the triangular Cholesky factor U or L.

      stptri(UPLO, 'Non-unit', N, AP, INFO );
      IF( INFO.GT.0 ) RETURN

      if ( UPPER ) {

         // Compute the product inv(U) * inv(U)**T.

         JJ = 0
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1
            JJ = JJ + J
            IF( J.GT.1 ) CALL SSPR( 'Upper', J-1, ONE, AP( JC ), 1, AP )
            AJJ = AP( JJ )
            sscal(J, AJJ, AP( JC ), 1 );
   10    CONTINUE

      } else {

         // Compute the product inv(L)**T * inv(L).

         JJ = 1
         for (J = 1; J <= N; J++) { // 20
            JJN = JJ + N - J + 1
            AP( JJ ) = SDOT( N-J+1, AP( JJ ), 1, AP( JJ ), 1 )
            IF( J.LT.N ) CALL STPMV( 'Lower', 'Transpose', 'Non-unit', N-J, AP( JJN ), AP( JJ+1 ), 1 )
            JJ = JJN
   20    CONTINUE
      }

      RETURN

      // End of SPPTRI

      }

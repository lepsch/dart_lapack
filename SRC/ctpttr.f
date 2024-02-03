      SUBROUTINE CTPTTR( UPLO, N, AP, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      int                I, J, K;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('CTPTTR', -INFO );
         RETURN
      }

      if ( LOWER ) {
         K = 0
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      } else {
         K = 0
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      }


      RETURN

      // End of CTPTTR

      }

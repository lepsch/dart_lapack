      SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AP( * )
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
         CALL XERBLA( 'ZTPTTR', -INFO )
         RETURN
      }

      if ( LOWER ) {
         K = 0
         DO J = 1, N
            DO I = J, N
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      } else {
         K = 0
         DO J = 1, N
            DO I = 1, J
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      }


      RETURN

      // End of ZTPTTR

      }

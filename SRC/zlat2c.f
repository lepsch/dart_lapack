      SUBROUTINE ZLAT2C( UPLO, N, A, LDA, SA, LDSA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDSA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            SA( LDSA, * )
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      bool               UPPER;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, CMPLX
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      bool               LSAME;
      // EXTERNAL SLAMCH, LSAME
      // ..
      // .. Executable Statements ..

      RMAX = SLAMCH( 'O' )
      UPPER = LSAME( UPLO, 'U' )
      if ( UPPER ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( ( DBLE( A( I, J ) ).LT.-RMAX ) .OR. ( DBLE( A( I, J ) ).GT.RMAX ) .OR. ( DIMAG( A( I, J ) ).LT.-RMAX ) .OR. ( DIMAG( A( I, J ) ).GT.RMAX ) ) {
                  INFO = 1
                  GO TO 50
               }
               SA( I, J ) = CMPLX( A( I, J ) )
   10       CONTINUE
   20    CONTINUE
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               if ( ( DBLE( A( I, J ) ).LT.-RMAX ) .OR. ( DBLE( A( I, J ) ).GT.RMAX ) .OR. ( DIMAG( A( I, J ) ).LT.-RMAX ) .OR. ( DIMAG( A( I, J ) ).GT.RMAX ) ) {
                  INFO = 1
                  GO TO 50
               }
               SA( I, J ) = CMPLX( A( I, J ) )
   30       CONTINUE
   40    CONTINUE
      }
   50 CONTINUE

      RETURN

      // End of ZLAT2C

      }

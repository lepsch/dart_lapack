      SUBROUTINE ZLAG2C( M, N, A, LDA, SA, LDSA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDSA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            SA( LDSA, * )
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, CMPLX
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Executable Statements ..

      RMAX = SLAMCH( 'O' )
      DO 20 J = 1, N
         DO 10 I = 1, M
            if ( ( DBLE( A( I, J ) ).LT.-RMAX ) .OR. ( DBLE( A( I, J ) ).GT.RMAX ) .OR. ( DIMAG( A( I, J ) ).LT.-RMAX ) .OR. ( DIMAG( A( I, J ) ).GT.RMAX ) ) {
               INFO = 1
               GO TO 30
            }
            SA( I, J ) = CMPLX( A( I, J ) )
   10    CONTINUE
   20 CONTINUE
      INFO = 0
   30 CONTINUE
      RETURN

      // End of ZLAG2C

      }

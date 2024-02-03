      SUBROUTINE ZLAT2C( UPLO, N, A, LDA, SA, LDSA, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDSA, N;
*     ..
*     .. Array Arguments ..
      COMPLEX            SA( LDSA, * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                I, J;
      double             RMAX;
      bool               UPPER;
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, CMPLX
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      bool               LSAME;
      EXTERNAL           SLAMCH, LSAME
*     ..
*     .. Executable Statements ..
*
      RMAX = SLAMCH( 'O' )
      UPPER = LSAME( UPLO, 'U' )
      IF( UPPER ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, J
               IF( ( DBLE( A( I, J ) ).LT.-RMAX ) .OR. ( DBLE( A( I, J ) ).GT.RMAX ) .OR. ( DIMAG( A( I, J ) ).LT.-RMAX ) .OR. ( DIMAG( A( I, J ) ).GT.RMAX ) ) THEN
                  INFO = 1
                  GO TO 50
               END IF
               SA( I, J ) = CMPLX( A( I, J ) )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = J, N
               IF( ( DBLE( A( I, J ) ).LT.-RMAX ) .OR. ( DBLE( A( I, J ) ).GT.RMAX ) .OR. ( DIMAG( A( I, J ) ).LT.-RMAX ) .OR. ( DIMAG( A( I, J ) ).GT.RMAX ) ) THEN
                  INFO = 1
                  GO TO 50
               END IF
               SA( I, J ) = CMPLX( A( I, J ) )
   30       CONTINUE
   40    CONTINUE
      END IF
   50 CONTINUE
*
      RETURN
*
*     End of ZLAT2C
*
      END

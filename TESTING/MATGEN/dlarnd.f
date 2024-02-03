      double           FUNCTION DLARND( IDIST, ISEED );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IDIST;
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE, TWO;
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      double             TWOPI;
      PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839D+0 )
*     ..
*     .. Local Scalars ..
      double             T1, T2;
*     ..
*     .. External Functions ..
      double             DLARAN;
      EXTERNAL           DLARAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, SQRT
*     ..
*     .. Executable Statements ..
*
*     Generate a real random number from a uniform (0,1) distribution
*
      T1 = DLARAN( ISEED )
*
      IF( IDIST.EQ.1 ) THEN
*
*        uniform (0,1)
*
         DLARND = T1
      ELSE IF( IDIST.EQ.2 ) THEN
*
*        uniform (-1,1)
*
         DLARND = TWO*T1 - ONE
      ELSE IF( IDIST.EQ.3 ) THEN
*
*        normal (0,1)
*
         T2 = DLARAN( ISEED )
         DLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 )
      END IF
      RETURN
*
*     End of DLARND
*
      END

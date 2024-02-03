      SUBROUTINE ZLARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX*16         X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      DOUBLE PRECISION   TWOPI
      PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IV
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   U( LV )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, EXP, LOG, MIN, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARUV
*     ..
*     .. Executable Statements ..
*
      DO 60 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
*
*        Call DLARUV to generate 2*IL real numbers from a uniform (0,1)
*        distribution (2*IL <= LV)
*
         CALL DLARUV( ISEED, 2*IL, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           Copy generated numbers
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = DCMPLX( U( 2*I-1 ), U( 2*I ) )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           Convert generated numbers to uniform (-1,1) distribution
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = DCMPLX( TWO*U( 2*I-1 )-ONE, TWO*U( 2*I )-ONE )
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           Convert generated numbers to normal (0,1) distribution
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )* EXP( DCMPLX( ZERO, TWOPI*U( 2*I ) ) )
   30       CONTINUE
         ELSE IF( IDIST.EQ.4 ) THEN
*
*           Convert generated numbers to complex numbers uniformly
*           distributed on the unit disk
*
            DO 40 I = 1, IL
               X( IV+I-1 ) = SQRT( U( 2*I-1 ) )* EXP( DCMPLX( ZERO, TWOPI*U( 2*I ) ) )
   40       CONTINUE
         ELSE IF( IDIST.EQ.5 ) THEN
*
*           Convert generated numbers to complex numbers uniformly
*           distributed on the unit circle
*
            DO 50 I = 1, IL
               X( IV+I-1 ) = EXP( DCMPLX( ZERO, TWOPI*U( 2*I ) ) )
   50       CONTINUE
         END IF
   60 CONTINUE
      RETURN
*
*     End of ZLARNV
*
      END

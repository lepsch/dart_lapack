      SUBROUTINE SGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      REAL               AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( * ), R( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               BIGNUM, RCMAX, RCMIN, SMLNUM
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
*     Get machine constants.
*
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
*     Compute row scale factors.
*
      DO 10 I = 1, M
         R( I ) = ZERO
   10 CONTINUE
*
*     Find the maximum element in each row.
*
      DO 30 J = 1, N
         DO 20 I = 1, M
            R( I ) = MAX( R( I ), ABS( A( I, J ) ) )
   20    CONTINUE
   30 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 50 I = 1, M
            IF( R( I ).EQ.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   50    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
*
*        Compute ROWCND = min(R(I)) / max(R(I))
*
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
*     Compute column scale factors
*
      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE
*
*     Find the maximum element in each column,
*     assuming the row scaling computed above.
*
      DO 90 J = 1, N
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), ABS( A( I, J ) )*R( I ) )
   80    CONTINUE
   90 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 110 J = 1, N
            IF( C( J ).EQ.ZERO ) THEN
               INFO = M + J
               RETURN
            END IF
  110    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
*
*        Compute COLCND = min(C(J)) / max(C(J))
*
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
      RETURN
*
*     End of SGEEQU
*
      END

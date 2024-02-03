      SUBROUTINE SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, N
*     ..
*     .. Array Arguments ..
      REAL               Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ALPHA, REALONE, REALZERO
      PARAMETER          ( ALPHA = 0.83E0, REALONE = 1.0E0, REALZERO = 0.0E0 )
      REAL               NEGONE, ONE, ZERO
      PARAMETER          ( NEGONE = -1.0E0, ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX
      REAL               EPS, NORM, NORM_NEW, SCL, SSQ
*     ..
*     .. External Functions ..
      REAL               SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SLASSQ, XERBLA
*     ..
*     .. Intrinsic Function ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test input arguments
*
      INFO = 0
      IF( M1 .LT. 0 ) THEN
         INFO = -1
      ELSE IF( M2 .LT. 0 ) THEN
         INFO = -2
      ELSE IF( N .LT. 0 ) THEN
         INFO = -3
      ELSE IF( INCX1 .LT. 1 ) THEN
         INFO = -5
      ELSE IF( INCX2 .LT. 1 ) THEN
         INFO = -7
      ELSE IF( LDQ1 .LT. MAX( 1, M1 ) ) THEN
         INFO = -9
      ELSE IF( LDQ2 .LT. MAX( 1, M2 ) ) THEN
         INFO = -11
      ELSE IF( LWORK .LT. N ) THEN
         INFO = -13
      END IF
*
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SORBDB6', -INFO )
         RETURN
      END IF
*
      EPS = SLAMCH( 'Precision' )
*
*     Compute the Euclidean norm of X
*
      SCL = REALZERO
      SSQ = REALZERO
      CALL SLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL SLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM = SCL * SQRT( SSQ )
*
*     First, project X onto the orthogonal complement of Q's column
*     space
*
      IF( M1 .EQ. 0 ) THEN
         DO I = 1, N
            WORK(I) = ZERO
         END DO
      ELSE
         CALL SGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 )
      END IF
*
      CALL SGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )
*
      CALL SGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL SGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 )
*
      SCL = REALZERO
      SSQ = REALZERO
      CALL SLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL SLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM_NEW = SCL * SQRT(SSQ)
*
*     If projection is sufficiently large in norm, then stop.
*     If projection is zero, then stop.
*     Otherwise, project again.
*
      IF( NORM_NEW .GE. ALPHA * NORM ) THEN
         RETURN
      END IF
*
      IF( NORM_NEW .LE. N * EPS * NORM ) THEN
         DO IX = 1, 1 + (M1-1)*INCX1, INCX1
           X1( IX ) = ZERO
         END DO
         DO IX = 1, 1 + (M2-1)*INCX2, INCX2
           X2( IX ) = ZERO
         END DO
         RETURN
      END IF
*
      NORM = NORM_NEW
*
      DO I = 1, N
         WORK(I) = ZERO
      END DO
*
      IF( M1 .EQ. 0 ) THEN
         DO I = 1, N
            WORK(I) = ZERO
         END DO
      ELSE
         CALL SGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 )
      END IF
*
      CALL SGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )
*
      CALL SGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL SGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 )
*
      SCL = REALZERO
      SSQ = REALZERO
      CALL SLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL SLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM_NEW = SCL * SQRT(SSQ)
*
*     If second projection is sufficiently large in norm, then do
*     nothing more. Alternatively, if it shrunk significantly, then
*     truncate it to zero.
*
      IF( NORM_NEW .LT. ALPHA * NORM ) THEN
         DO IX = 1, 1 + (M1-1)*INCX1, INCX1
            X1(IX) = ZERO
         END DO
         DO IX = 1, 1 + (M2-1)*INCX2, INCX2
            X2(IX) = ZERO
         END DO
      END IF
*
      RETURN
*
*     End of SORBDB6
*
      END

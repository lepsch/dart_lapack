      SUBROUTINE DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, N;
      // ..
      // .. Array Arguments ..
      double             Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ALPHA, REALONE, REALZERO;
      PARAMETER          ( ALPHA = 0.83D0, REALONE = 1.0D0, REALZERO = 0.0D0 )
      double             NEGONE, ONE, ZERO;
      PARAMETER          ( NEGONE = -1.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      int                I, IX;
      double             EPS, NORM, NORM_NEW, SCL, SSQ;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DLASSQ, XERBLA
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test input arguments

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

      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DORBDB6', -INFO )
         RETURN
      END IF

      EPS = DLAMCH( 'Precision' )

      // Compute the Euclidean norm of X

      SCL = REALZERO
      SSQ = REALZERO
      CALL DLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL DLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM = SCL * SQRT( SSQ )

      // First, project X onto the orthogonal complement of Q's column
      // space

      IF( M1 .EQ. 0 ) THEN
         DO I = 1, N
            WORK(I) = ZERO
         END DO
      ELSE
         CALL DGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 )
      END IF

      CALL DGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )

      CALL DGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL DGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 )

      SCL = REALZERO
      SSQ = REALZERO
      CALL DLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL DLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM_NEW = SCL * SQRT(SSQ)

      // If projection is sufficiently large in norm, then stop.
      // If projection is zero, then stop.
      // Otherwise, project again.

      IF( NORM_NEW .GE. ALPHA * NORM ) THEN
         RETURN
      END IF

      IF( NORM_NEW .LE. N * EPS * NORM ) THEN
         DO IX = 1, 1 + (M1-1)*INCX1, INCX1
           X1( IX ) = ZERO
         END DO
         DO IX = 1, 1 + (M2-1)*INCX2, INCX2
           X2( IX ) = ZERO
         END DO
         RETURN
      END IF

      NORM = NORM_NEW

      DO I = 1, N
         WORK(I) = ZERO
      END DO

      IF( M1 .EQ. 0 ) THEN
         DO I = 1, N
            WORK(I) = ZERO
         END DO
      ELSE
         CALL DGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 )
      END IF

      CALL DGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )

      CALL DGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL DGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 )

      SCL = REALZERO
      SSQ = REALZERO
      CALL DLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL DLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM_NEW = SCL * SQRT(SSQ)

      // If second projection is sufficiently large in norm, then do
      // nothing more. Alternatively, if it shrunk significantly, then
     t // runcate it to zero.

      IF( NORM_NEW .LT. ALPHA * NORM ) THEN
         DO IX = 1, 1 + (M1-1)*INCX1, INCX1
            X1(IX) = ZERO
         END DO
         DO IX = 1, 1 + (M2-1)*INCX2, INCX2
            X2(IX) = ZERO
         END DO
      END IF

      RETURN

      // End of DORBDB6

      }

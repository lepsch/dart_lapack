      SUBROUTINE ZUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ALPHA, REALONE, REALZERO;
      const              ALPHA = 0.83D0, REALONE = 1.0D0, REALZERO = 0.0D0 ;
      COMPLEX*16         NEGONE, ONE, ZERO
      const              NEGONE = (-1.0D0,0.0D0), ONE = (1.0D0,0.0D0), ZERO = (0.0D0,0.0D0) ;
      // ..
      // .. Local Scalars ..
      int                I, IX;
      double             EPS, NORM, NORM_NEW, SCL, SSQ;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMV, ZLASSQ, XERBLA
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      if ( M1 .LT. 0 ) {
         INFO = -1
      } else if ( M2 .LT. 0 ) {
         INFO = -2
      } else if ( N .LT. 0 ) {
         INFO = -3
      } else if ( INCX1 .LT. 1 ) {
         INFO = -5
      } else if ( INCX2 .LT. 1 ) {
         INFO = -7
      } else if ( LDQ1 .LT. MAX( 1, M1 ) ) {
         INFO = -9
      } else if ( LDQ2 .LT. MAX( 1, M2 ) ) {
         INFO = -11
      } else if ( LWORK .LT. N ) {
         INFO = -13
      }

      if ( INFO .NE. 0 ) {
         CALL XERBLA( 'ZUNBDB6', -INFO )
         RETURN
      }

      EPS = DLAMCH( 'Precision' )

      // Compute the Euclidean norm of X

      SCL = REALZERO
      SSQ = REALZERO
      CALL ZLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL ZLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM = SCL * SQRT( SSQ )

      // First, project X onto the orthogonal complement of Q's column
      // space

      if ( M1 .EQ. 0 ) {
         DO I = 1, N
            WORK(I) = ZERO
         END DO
      } else {
         CALL ZGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 )
      }

      CALL ZGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )

      CALL ZGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL ZGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 )

      SCL = REALZERO
      SSQ = REALZERO
      CALL ZLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL ZLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM_NEW = SCL * SQRT(SSQ)

      // If projection is sufficiently large in norm, then stop.
      // If projection is zero, then stop.
      // Otherwise, project again.

      if ( NORM_NEW .GE. ALPHA * NORM ) {
         RETURN
      }

      if ( NORM_NEW .LE. N * EPS * NORM ) {
         DO IX = 1, 1 + (M1-1)*INCX1, INCX1
           X1( IX ) = ZERO
         END DO
         DO IX = 1, 1 + (M2-1)*INCX2, INCX2
           X2( IX ) = ZERO
         END DO
         RETURN
      }

      NORM = NORM_NEW

      DO I = 1, N
         WORK(I) = ZERO
      END DO

      if ( M1 .EQ. 0 ) {
         DO I = 1, N
            WORK(I) = ZERO
         END DO
      } else {
         CALL ZGEMV( 'C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 )
      }

      CALL ZGEMV( 'C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 )

      CALL ZGEMV( 'N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL ZGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 )

      SCL = REALZERO
      SSQ = REALZERO
      CALL ZLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL ZLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM_NEW = SCL * SQRT(SSQ)

      // If second projection is sufficiently large in norm, then do
      // nothing more. Alternatively, if it shrunk significantly, then
     t // runcate it to zero.

      if ( NORM_NEW .LT. ALPHA * NORM ) {
         DO IX = 1, 1 + (M1-1)*INCX1, INCX1
            X1(IX) = ZERO
         END DO
         DO IX = 1, 1 + (M2-1)*INCX2, INCX2
            X2(IX) = ZERO
         END DO
      }

      RETURN

      // End of ZUNBDB6

      }

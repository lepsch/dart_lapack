      SUBROUTINE DORBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, INFO )

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
      double             REALZERO;
      PARAMETER          ( REALZERO = 0.0D0 )
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      int                CHILDINFO, I, J;
      double             EPS, NORM, SCL, SSQ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ, DORBDB6, DSCAL, XERBLA
      // ..
      // .. External Functions ..
      double             DLAMCH, DNRM2;
      // EXTERNAL DLAMCH, DNRM2
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
         CALL XERBLA( 'DORBDB5', -INFO )
         RETURN
      END IF

      EPS = DLAMCH( 'Precision' )

      // Project X onto the orthogonal complement of Q if X is nonzero

      SCL = REALZERO
      SSQ = REALZERO
      CALL DLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL DLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM = SCL * SQRT( SSQ )

      IF( NORM .GT. N * EPS ) THEN
         // Scale vector to unit norm to avoid problems in the caller code.
         // Computing the reciprocal is undesirable but
          // * xLASCL cannot be used because of the vector increments and
          // * the round-off error has a negligible impact on
            // orthogonalization.
         CALL DSCAL( M1, ONE / NORM, X1, INCX1 )
         CALL DSCAL( M2, ONE / NORM, X2, INCX2 )
         CALL DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )

         // If the projection is nonzero, then return

         IF( DNRM2(M1,X1,INCX1) .NE. REALZERO .OR. DNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         END IF
      END IF

      // Project each standard basis vector e_1,...,e_M1 in turn, stopping
      // when a nonzero projection is found

      DO I = 1, M1
         DO J = 1, M1
            X1(J) = ZERO
         END DO
         X1(I) = ONE
         DO J = 1, M2
            X2(J) = ZERO
         END DO
         CALL DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )          IF( DNRM2(M1,X1,INCX1) .NE. REALZERO .OR. DNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         END IF
      END DO

      // Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
      // stopping when a nonzero projection is found

      DO I = 1, M2
         DO J = 1, M1
            X1(J) = ZERO
         END DO
         DO J = 1, M2
            X2(J) = ZERO
         END DO
         X2(I) = ONE
         CALL DORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )          IF( DNRM2(M1,X1,INCX1) .NE. REALZERO .OR. DNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         END IF
      END DO

      RETURN

      // End of DORBDB5

      END

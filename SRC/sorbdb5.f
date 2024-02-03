      SUBROUTINE SORBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, N;
      // ..
      // .. Array Arguments ..
      REAL               Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               REALZERO
      const              REALZERO = 0.0E0 ;
      REAL               ONE, ZERO
      const              ONE = 1.0E0, ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                CHILDINFO, I, J;
      REAL               EPS, NORM, SCL, SSQ
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ, SORBDB6, SSCAL, XERBLA
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SNRM2
      // EXTERNAL SLAMCH, SNRM2
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
         xerbla('SORBDB5', -INFO );
         RETURN
      }

      EPS = SLAMCH( 'Precision' )

      // Project X onto the orthogonal complement of Q if X is nonzero

      SCL = REALZERO
      SSQ = REALZERO
      slassq(M1, X1, INCX1, SCL, SSQ );
      slassq(M2, X2, INCX2, SCL, SSQ );
      NORM = SCL * SQRT( SSQ )

      if ( NORM .GT. N * EPS ) {
         // Scale vector to unit norm to avoid problems in the caller code.
         // Computing the reciprocal is undesirable but
          // * xLASCL cannot be used because of the vector increments and
          // * the round-off error has a negligible impact on
            // orthogonalization.
         sscal(M1, ONE / NORM, X1, INCX1 );
         sscal(M2, ONE / NORM, X2, INCX2 );
         sorbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO );

         // If the projection is nonzero, then return

         if ( SNRM2(M1,X1,INCX1) .NE. REALZERO .OR. SNRM2(M2,X2,INCX2) .NE. REALZERO ) {
            RETURN
         }
      }

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
         CALL SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )          IF( SNRM2(M1,X1,INCX1) .NE. REALZERO .OR. SNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         }
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
         CALL SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )          IF( SNRM2(M1,X1,INCX1) .NE. REALZERO .OR. SNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         }
      END DO

      RETURN

      // End of SORBDB5

      }

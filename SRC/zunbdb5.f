      SUBROUTINE ZUNBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
     $                    LDQ2, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2,
     $                   N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   REALZERO
      PARAMETER          ( REALZERO = 0.0D0 )
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = (1.0D0,0.0D0), ZERO = (0.0D0,0.0D0) )
*     ..
*     .. Local Scalars ..
      INTEGER            CHILDINFO, I, J
      DOUBLE PRECISION   EPS, NORM, SCL, SSQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ, ZUNBDB6, ZSCAL, XERBLA
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZNRM2
      EXTERNAL           DLAMCH, DZNRM2
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
         CALL XERBLA( 'ZUNBDB5', -INFO )
         RETURN
      END IF
*
      EPS = DLAMCH( 'Precision' )
*
*     Project X onto the orthogonal complement of Q if X is nonzero
*
      SCL = REALZERO
      SSQ = REALZERO
      CALL ZLASSQ( M1, X1, INCX1, SCL, SSQ )
      CALL ZLASSQ( M2, X2, INCX2, SCL, SSQ )
      NORM = SCL * SQRT( SSQ )
*
      IF( NORM .GT. N * EPS ) THEN
*        Scale vector to unit norm to avoid problems in the caller code.
*        Computing the reciprocal is undesirable but
*         * xLASCL cannot be used because of the vector increments and
*         * the round-off error has a negligible impact on
*           orthogonalization.
         CALL ZSCAL( M1, ONE / NORM, X1, INCX1 )
         CALL ZSCAL( M2, ONE / NORM, X2, INCX2 )
         CALL ZUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
     $              LDQ2, WORK, LWORK, CHILDINFO )
*
*        If the projection is nonzero, then return
*
         IF( DZNRM2(M1,X1,INCX1) .NE. REALZERO
     $       .OR. DZNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         END IF
      END IF
*
*     Project each standard basis vector e_1,...,e_M1 in turn, stopping
*     when a nonzero projection is found
*
      DO I = 1, M1
         DO J = 1, M1
            X1(J) = ZERO
         END DO
         X1(I) = ONE
         DO J = 1, M2
            X2(J) = ZERO
         END DO
         CALL ZUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
     $                 LDQ2, WORK, LWORK, CHILDINFO )
         IF( DZNRM2(M1,X1,INCX1) .NE. REALZERO
     $       .OR. DZNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         END IF
      END DO
*
*     Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
*     stopping when a nonzero projection is found
*
      DO I = 1, M2
         DO J = 1, M1
            X1(J) = ZERO
         END DO
         DO J = 1, M2
            X2(J) = ZERO
         END DO
         X2(I) = ONE
         CALL ZUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2,
     $                 LDQ2, WORK, LWORK, CHILDINFO )
         IF( DZNRM2(M1,X1,INCX1) .NE. REALZERO
     $       .OR. DZNRM2(M2,X2,INCX2) .NE. REALZERO ) THEN
            RETURN
         END IF
      END DO
*
      RETURN
*
*     End of ZUNBDB5
*
      END

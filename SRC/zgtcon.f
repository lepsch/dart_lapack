      SUBROUTINE ZGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         D( * ), DL( * ), DU( * ), DU2( * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      int                I, KASE, KASE1;
      double             AINVNM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGTTRS, ZLACN2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGTCON', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      // Check that D(1:N) is non-zero.
*
      DO 10 I = 1, N
         IF( D( I ).EQ.DCMPLX( ZERO ) ) RETURN
   10 CONTINUE
*
      AINVNM = ZERO
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   20 CONTINUE
      CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
            // Multiply by inv(U)*inv(L).
*
            CALL ZGTTRS( 'No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO )
         ELSE
*
            // Multiply by inv(L**H)*inv(U**H).
*
            CALL ZGTTRS( 'Conjugate transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO )
         END IF
         GO TO 20
      END IF
*
      // Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
*
      // End of ZGTCON
*
      END

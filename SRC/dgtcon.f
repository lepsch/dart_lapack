      SUBROUTINE DGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             D( * ), DL( * ), DU( * ), DU2( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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
      // EXTERNAL DGTTRS, DLACN2, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

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
         CALL XERBLA( 'DGTCON', -INFO )
         RETURN
      END IF

      // Quick return if possible

      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF

      // Check that D(1:N) is non-zero.

      DO 10 I = 1, N
         IF( D( I ).EQ.ZERO ) RETURN
   10 CONTINUE

      AINVNM = ZERO
      IF( ONENRM ) THEN
         KASE1 = 1
      } else {
         KASE1 = 2
      END IF
      KASE = 0
   20 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN

            // Multiply by inv(U)*inv(L).

            CALL DGTTRS( 'No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO )
         } else {

            // Multiply by inv(L**T)*inv(U**T).

            CALL DGTTRS( 'Transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO )
         END IF
         GO TO 20
      END IF

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      RETURN

      // End of DGTCON

      }

      SUBROUTINE ZSYCON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), E( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, KASE;
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
      // EXTERNAL ZLACN2, ZSYTRS_3, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYCON_3', -INFO )
         RETURN
      END IF

      // Quick return if possible

      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.LE.ZERO ) THEN
         RETURN
      END IF

      // Check that the diagonal matrix D is nonsingular.

      IF( UPPER ) THEN

         // Upper triangular storage: examine D from bottom to top

         DO I = N, 1, -1
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.CZERO ) RETURN
         END DO
      ELSE

         // Lower triangular storage: examine D from top to bottom.

         DO I = 1, N
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.CZERO ) RETURN
         END DO
      END IF

      // Estimate the 1-norm of the inverse.

      KASE = 0
   30 CONTINUE
      CALL ZLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN

         // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

         CALL ZSYTRS_3( UPLO, N, 1, A, LDA, E, IPIV, WORK, N, INFO )
         GO TO 30
      END IF

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      RETURN

      // End of ZSYCON_3

      }

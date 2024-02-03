      SUBROUTINE SGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               D( * ), DL( * ), DU( * ), DU2( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      int                I, KASE, KASE1;
      REAL               AINVNM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGTTRS, SLACN2, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      if ( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( ANORM.LT.ZERO ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('SGTCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N.EQ.0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM.EQ.ZERO ) {
         RETURN
      }

      // Check that D(1:N) is non-zero.

      DO 10 I = 1, N
         IF( D( I ).EQ.ZERO ) RETURN
   10 CONTINUE

      AINVNM = ZERO
      if ( ONENRM ) {
         KASE1 = 1
      } else {
         KASE1 = 2
      }
      KASE = 0
   20 CONTINUE
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.KASE1 ) {

            // Multiply by inv(U)*inv(L).

            sgttrs('No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO );
         } else {

            // Multiply by inv(L**T)*inv(U**T).

            sgttrs('Transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO );
         }
         GO TO 20
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      RETURN

      // End of SGTCON

      }

      SUBROUTINE ZPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( 0: * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, NORMALTRANSR;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTFSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPFTRS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      // start execution: there are two triangular solves

      IF( LOWER ) THEN
         CALL ZTFSM( TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B, LDB )          CALL ZTFSM( TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B, LDB )
      } else {
         CALL ZTFSM( TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B, LDB )          CALL ZTFSM( TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B, LDB )
      END IF

      RETURN

      // End of ZPFTRS

      }

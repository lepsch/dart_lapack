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
      if ( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) {
         INFO = -1
      } else if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         xerbla('ZPFTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0 .OR. NRHS.EQ.0) RETURN;

      // start execution: there are two triangular solves

      if ( LOWER ) {
         ztfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B, LDB )          CALL ZTFSM( TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B, LDB );
      } else {
         ztfsm(TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B, LDB )          CALL ZTFSM( TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B, LDB );
      }

      RETURN

      // End of ZPFTRS

      }

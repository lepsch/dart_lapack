      SUBROUTINE CTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            AP( * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JC;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( NRHS.LT.0 ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('CTPTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      // Check for singularity.

      if ( NOUNIT ) {
         if ( UPPER ) {
            JC = 1
            for (INFO = 1; INFO <= N; INFO++) { // 10
               IF( AP( JC+INFO-1 ).EQ.ZERO ) RETURN
               JC = JC + INFO
            } // 10
         } else {
            JC = 1
            for (INFO = 1; INFO <= N; INFO++) { // 20
               IF( AP( JC ).EQ.ZERO ) RETURN
               JC = JC + N - INFO + 1
            } // 20
         }
      }
      INFO = 0

      // Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.

      for (J = 1; J <= NRHS; J++) { // 30
         ctpsv(UPLO, TRANS, DIAG, N, AP, B( 1, J ), 1 );
      } // 30

      RETURN

      // End of CTPTRS

      }

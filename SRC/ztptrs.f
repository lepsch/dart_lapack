      SUBROUTINE ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
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
      // EXTERNAL XERBLA, ZTPSV
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
         xerbla('ZTPTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Check for singularity.

      if ( NOUNIT ) {
         if ( UPPER ) {
            JC = 1
            DO 10 INFO = 1, N
               IF( AP( JC+INFO-1 ).EQ.ZERO ) RETURN
               JC = JC + INFO
   10       CONTINUE
         } else {
            JC = 1
            DO 20 INFO = 1, N
               IF( AP( JC ).EQ.ZERO ) RETURN
               JC = JC + N - INFO + 1
   20       CONTINUE
         }
      }
      INFO = 0

      // Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.

      DO 30 J = 1, NRHS
         ztpsv(UPLO, TRANS, DIAG, N, AP, B( 1, J ), 1 );
   30 CONTINUE

      RETURN

      // End of ZTPTRS

      }

      SUBROUTINE ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRAN;
      int                ITRANS, J, JB, NB;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGTTS2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOTRAN = ( TRANS.EQ.'N' .OR. TRANS.EQ.'n' )
      if ( .NOT.NOTRAN .AND. .NOT.( TRANS.EQ.'T' .OR. TRANS.EQ. 't' ) .AND. .NOT.( TRANS.EQ.'C' .OR. TRANS.EQ.'c' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDB.LT.MAX( N, 1 ) ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGTTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      // Decode TRANS

      if ( NOTRAN ) {
         ITRANS = 0
      } else if ( TRANS.EQ.'T' .OR. TRANS.EQ.'t' ) {
         ITRANS = 1
      } else {
         ITRANS = 2
      }

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS.EQ.1 ) {
         NB = 1
      } else {
         NB = MAX( 1, ILAENV( 1, 'ZGTTRS', TRANS, N, NRHS, -1, -1 ) )
      }

      if ( NB.GE.NRHS ) {
         zgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB );
      } else {
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            zgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), LDB );
         } // 10
      }

      // End of ZGTTRS

      }

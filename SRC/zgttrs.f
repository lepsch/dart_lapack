      SUBROUTINE ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               NOTRAN;
      int                ITRANS, J, JB, NB
*     ..
*     .. External Functions ..
      int                ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGTTS2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOTRAN = ( TRANS.EQ.'N' .OR. TRANS.EQ.'n' )
      IF( .NOT.NOTRAN .AND. .NOT.( TRANS.EQ.'T' .OR. TRANS.EQ. 't' ) .AND. .NOT.( TRANS.EQ.'C' .OR. TRANS.EQ.'c' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
*     Decode TRANS
*
      IF( NOTRAN ) THEN
         ITRANS = 0
      ELSE IF( TRANS.EQ.'T' .OR. TRANS.EQ.'t' ) THEN
         ITRANS = 1
      ELSE
         ITRANS = 2
      END IF
*
*     Determine the number of right-hand sides to solve at a time.
*
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'ZGTTRS', TRANS, N, NRHS, -1, -1 ) )
      END IF
*
      IF( NB.GE.NRHS ) THEN
         CALL ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL ZGTTS2( ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), LDB )
   10    CONTINUE
      END IF
*
*     End of ZGTTRS
*
      END

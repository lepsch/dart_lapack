      SUBROUTINE CTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, N, NRHS;
*     ..
*     .. Array Arguments ..
      COMPLEX            AP( * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JC;
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTPSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Check for singularity.
*
      IF( NOUNIT ) THEN
         IF( UPPER ) THEN
            JC = 1
            DO 10 INFO = 1, N
               IF( AP( JC+INFO-1 ).EQ.ZERO ) RETURN
               JC = JC + INFO
   10       CONTINUE
         ELSE
            JC = 1
            DO 20 INFO = 1, N
               IF( AP( JC ).EQ.ZERO ) RETURN
               JC = JC + N - INFO + 1
   20       CONTINUE
         END IF
      END IF
      INFO = 0
*
*     Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.
*
      DO 30 J = 1, NRHS
         CALL CTPSV( UPLO, TRANS, DIAG, N, AP, B( 1, J ), 1 )
   30 CONTINUE
*
      RETURN
*
*     End of CTPTRS
*
      END

      SUBROUTINE ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
      COMPLEX*16         B( LDB, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               UPPER;
      int                IUPLO, J, JB, NB
*     ..
*     .. External Functions ..
      int                ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZPTTS2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      UPPER = ( UPLO.EQ.'U' .OR. UPLO.EQ.'u' )
      IF( .NOT.UPPER .AND. .NOT.( UPLO.EQ.'L' .OR. UPLO.EQ.'l' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
*     Determine the number of right-hand sides to solve at a time.
*
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'ZPTTRS', UPLO, N, NRHS, -1, -1 ) )
      END IF
*
*     Decode UPLO
*
      IF( UPPER ) THEN
         IUPLO = 1
      ELSE
         IUPLO = 0
      END IF
*
      IF( NB.GE.NRHS ) THEN
         CALL ZPTTS2( IUPLO, N, NRHS, D, E, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL ZPTTS2( IUPLO, N, JB, D, E, B( 1, J ), LDB )
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of ZPTTRS
*
      END

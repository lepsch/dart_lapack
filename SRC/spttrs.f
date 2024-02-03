      SUBROUTINE SPTTRS( N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                J, JB, NB;
*     ..
*     .. External Functions ..
      int                ILAENV;
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           SPTTS2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPTTRS', -INFO )
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
         NB = MAX( 1, ILAENV( 1, 'SPTTRS', ' ', N, NRHS, -1, -1 ) )
      END IF
*
      IF( NB.GE.NRHS ) THEN
         CALL SPTTS2( N, NRHS, D, E, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL SPTTS2( N, JB, D, E, B( 1, J ), LDB )
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of SPTTRS
*
      END

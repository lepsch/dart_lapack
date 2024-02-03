      SUBROUTINE SPPTRI( UPLO, N, AP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INFO, N
*     ..
*     .. Array Arguments ..
      REAL               AP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      int                J, JC, JJ, JJN
      REAL               AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSPR, STPMV, STPTRI, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPPTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Invert the triangular Cholesky factor U or L.
*
      CALL STPTRI( UPLO, 'Non-unit', N, AP, INFO )
      IF( INFO.GT.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Compute the product inv(U) * inv(U)**T.
*
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
            IF( J.GT.1 ) CALL SSPR( 'Upper', J-1, ONE, AP( JC ), 1, AP )
            AJJ = AP( JJ )
            CALL SSCAL( J, AJJ, AP( JC ), 1 )
   10    CONTINUE
*
      ELSE
*
*        Compute the product inv(L)**T * inv(L).
*
         JJ = 1
         DO 20 J = 1, N
            JJN = JJ + N - J + 1
            AP( JJ ) = SDOT( N-J+1, AP( JJ ), 1, AP( JJ ), 1 )
            IF( J.LT.N ) CALL STPMV( 'Lower', 'Transpose', 'Non-unit', N-J, AP( JJN ), AP( JJ+1 ), 1 )
            JJ = JJN
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of SPPTRI
*
      END

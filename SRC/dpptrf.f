      SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
*     ..
*     .. Array Arguments ..
      double             AP( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ;
      double             AJJ;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DDOT;
      EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSPR, DTPSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC SQRT
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
         CALL XERBLA( 'DPPTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U**T*U.
*
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
*
*           Compute elements 1:J-1 of column J.
*
            IF( J.GT.1 ) CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', J-1, AP, AP( JC ), 1 )
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = AP( JJ ) - DDOT( J-1, AP( JC ), 1, AP( JC ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AP( JJ ) = SQRT( AJJ )
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L**T.
*
         JJ = 1
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = AP( JJ )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AP( JJ ) = AJJ
*
*           Compute elements J+1:N of column J and update the trailing
*           submatrix.
*
            IF( J.LT.N ) THEN
               CALL DSCAL( N-J, ONE / AJJ, AP( JJ+1 ), 1 )
               CALL DSPR( 'Lower', N-J, -ONE, AP( JJ+1 ), 1, AP( JJ+N-J+1 ) )
               JJ = JJ + N - J + 1
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of DPPTRF
*
      END

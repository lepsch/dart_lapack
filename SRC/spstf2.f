      SUBROUTINE SPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               TOL
      int                INFO, LDA, N, RANK
      String             UPLO;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( 2*N )
      int                PIV( N )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      REAL               AJJ, SSTOP, STEMP
      int                I, ITEMP, J, PVT
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      LOGICAL            LSAME, SISNAN
      EXTERNAL           SLAMCH, LSAME, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SSCAL, SSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT, MAXLOC
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPSTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Initialize PIV
*
      DO 100 I = 1, N
         PIV( I ) = I
  100 CONTINUE
*
*     Compute stopping value
*
      PVT = 1
      AJJ = A( PVT, PVT )
      DO I = 2, N
         IF( A( I, I ).GT.AJJ ) THEN
            PVT = I
            AJJ = A( PVT, PVT )
         END IF
      END DO
      IF( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) THEN
         RANK = 0
         INFO = 1
         GO TO 170
      END IF
*
*     Compute stopping value if not supplied
*
      IF( TOL.LT.ZERO ) THEN
         SSTOP = N * SLAMCH( 'Epsilon' ) * AJJ
      ELSE
         SSTOP = TOL
      END IF
*
*     Set first half of WORK to zero, holds dot products
*
      DO 110 I = 1, N
         WORK( I ) = 0
  110 CONTINUE
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization P**T * A * P = U**T * U
*
         DO 130 J = 1, N
*
*        Find pivot, test for exit, else swap rows and columns
*        Update dot products, compute possible pivots which are
*        stored in the second half of WORK
*
            DO 120 I = J, N
*
               IF( J.GT.1 ) THEN
                  WORK( I ) = WORK( I ) + A( J-1, I )**2
               END IF
               WORK( N+I ) = A( I, I ) - WORK( I )
*
  120       CONTINUE
*
            IF( J.GT.1 ) THEN
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.SSTOP.OR.SISNAN( AJJ ) ) THEN
                  A( J, J ) = AJJ
                  GO TO 160
               END IF
            END IF
*
            IF( J.NE.PVT ) THEN
*
*              Pivot OK, so can now swap pivot rows and columns
*
               A( PVT, PVT ) = A( J, J )
               CALL SSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
               IF( PVT.LT.N ) CALL SSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA )
               CALL SSWAP( PVT-J-1, A( J, J+1 ), LDA, A( J+1, PVT ), 1 )
*
*              Swap dot products and PIV
*
               STEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = STEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF
*
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J
*
            IF( J.LT.N ) THEN
               CALL SGEMV( 'Trans', J-1, N-J, -ONE, A( 1, J+1 ), LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL SSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
*
  130    CONTINUE
*
      ELSE
*
*        Compute the Cholesky factorization P**T * A * P = L * L**T
*
         DO 150 J = 1, N
*
*        Find pivot, test for exit, else swap rows and columns
*        Update dot products, compute possible pivots which are
*        stored in the second half of WORK
*
            DO 140 I = J, N
*
               IF( J.GT.1 ) THEN
                  WORK( I ) = WORK( I ) + A( I, J-1 )**2
               END IF
               WORK( N+I ) = A( I, I ) - WORK( I )
*
  140       CONTINUE
*
            IF( J.GT.1 ) THEN
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.SSTOP.OR.SISNAN( AJJ ) ) THEN
                  A( J, J ) = AJJ
                  GO TO 160
               END IF
            END IF
*
            IF( J.NE.PVT ) THEN
*
*              Pivot OK, so can now swap pivot rows and columns
*
               A( PVT, PVT ) = A( J, J )
               CALL SSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
               IF( PVT.LT.N ) CALL SSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
               CALL SSWAP( PVT-J-1, A( J+1, J ), 1, A( PVT, J+1 ), LDA )
*
*              Swap dot products and PIV
*
               STEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = STEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF
*
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J
*
            IF( J.LT.N ) THEN
               CALL SGEMV( 'No Trans', N-J, J-1, -ONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL SSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
*
  150    CONTINUE
*
      END IF
*
*     Ran to completion, A has full rank
*
      RANK = N
*
      GO TO 170
  160 CONTINUE
*
*     Rank is number of steps completed.  Set INFO = 1 to signal
*     that the factorization cannot be used to solve a system.
*
      RANK = J - 1
      INFO = 1
*
  170 CONTINUE
      RETURN
*
*     End of SPSTF2
*
      END

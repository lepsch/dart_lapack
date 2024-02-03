      SUBROUTINE DPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   TOL
      int                INFO, LDA, N, RANK
      CHARACTER          UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( 2*N )
      int                PIV( N )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AJJ, DSTOP, DTEMP
      int                I, ITEMP, J, JB, K, NB, PVT
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      int                ILAENV
      LOGICAL            LSAME, DISNAN
      EXTERNAL           DLAMCH, ILAENV, LSAME, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DPSTF2, DSCAL, DSWAP, DSYRK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT, MAXLOC
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
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPSTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Get block size
*
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DPSTF2( UPLO, N, A( 1, 1 ), LDA, PIV, RANK, TOL, WORK, INFO )
         GO TO 200
*
      ELSE
*
*     Initialize PIV
*
         DO 100 I = 1, N
            PIV( I ) = I
  100    CONTINUE
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
         IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
            RANK = 0
            INFO = 1
            GO TO 200
         END IF
*
*     Compute stopping value if not supplied
*
         IF( TOL.LT.ZERO ) THEN
            DSTOP = N * DLAMCH( 'Epsilon' ) * AJJ
         ELSE
            DSTOP = TOL
         END IF
*
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization P**T * A * P = U**T * U
*
            DO 140 K = 1, N, NB
*
*              Account for last block not being NB wide
*
               JB = MIN( NB, N-K+1 )
*
*              Set relevant part of first half of WORK to zero,
*              holds dot products
*
               DO 110 I = K, N
                  WORK( I ) = 0
  110          CONTINUE
*
               DO 130 J = K, K + JB - 1
*
*              Find pivot, test for exit, else swap rows and columns
*              Update dot products, compute possible pivots which are
*              stored in the second half of WORK
*
                  DO 120 I = J, N
*
                     IF( J.GT.K ) THEN
                        WORK( I ) = WORK( I ) + A( J-1, I )**2
                     END IF
                     WORK( N+I ) = A( I, I ) - WORK( I )
*
  120             CONTINUE
*
                  IF( J.GT.1 ) THEN
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                     PVT = ITEMP + J - 1
                     AJJ = WORK( N+PVT )
                     IF( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) THEN
                        A( J, J ) = AJJ
                        GO TO 190
                     END IF
                  END IF
*
                  IF( J.NE.PVT ) THEN
*
*                    Pivot OK, so can now swap pivot rows and columns
*
                     A( PVT, PVT ) = A( J, J )
                     CALL DSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
                     IF( PVT.LT.N ) CALL DSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA )
                     CALL DSWAP( PVT-J-1, A( J, J+1 ), LDA, A( J+1, PVT ), 1 )
*
*                    Swap dot products and PIV
*
                     DTEMP = WORK( J )
                     WORK( J ) = WORK( PVT )
                     WORK( PVT ) = DTEMP
                     ITEMP = PIV( PVT )
                     PIV( PVT ) = PIV( J )
                     PIV( J ) = ITEMP
                  END IF
*
                  AJJ = SQRT( AJJ )
                  A( J, J ) = AJJ
*
*                 Compute elements J+1:N of row J.
*
                  IF( J.LT.N ) THEN
                     CALL DGEMV( 'Trans', J-K, N-J, -ONE, A( K, J+1 ), LDA, A( K, J ), 1, ONE, A( J, J+1 ), LDA )
                     CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
                  END IF
*
  130          CONTINUE
*
*              Update trailing matrix, J already incremented
*
               IF( K+JB.LE.N ) THEN
                  CALL DSYRK( 'Upper', 'Trans', N-J+1, JB, -ONE, A( K, J ), LDA, ONE, A( J, J ), LDA )
               END IF
*
  140       CONTINUE
*
         ELSE
*
*        Compute the Cholesky factorization P**T * A * P = L * L**T
*
            DO 180 K = 1, N, NB
*
*              Account for last block not being NB wide
*
               JB = MIN( NB, N-K+1 )
*
*              Set relevant part of first half of WORK to zero,
*              holds dot products
*
               DO 150 I = K, N
                  WORK( I ) = 0
  150          CONTINUE
*
               DO 170 J = K, K + JB - 1
*
*              Find pivot, test for exit, else swap rows and columns
*              Update dot products, compute possible pivots which are
*              stored in the second half of WORK
*
                  DO 160 I = J, N
*
                     IF( J.GT.K ) THEN
                        WORK( I ) = WORK( I ) + A( I, J-1 )**2
                     END IF
                     WORK( N+I ) = A( I, I ) - WORK( I )
*
  160             CONTINUE
*
                  IF( J.GT.1 ) THEN
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                     PVT = ITEMP + J - 1
                     AJJ = WORK( N+PVT )
                     IF( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) THEN
                        A( J, J ) = AJJ
                        GO TO 190
                     END IF
                  END IF
*
                  IF( J.NE.PVT ) THEN
*
*                    Pivot OK, so can now swap pivot rows and columns
*
                     A( PVT, PVT ) = A( J, J )
                     CALL DSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
                     IF( PVT.LT.N ) CALL DSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
                     CALL DSWAP( PVT-J-1, A( J+1, J ), 1, A( PVT, J+1 ), LDA )
*
*                    Swap dot products and PIV
*
                     DTEMP = WORK( J )
                     WORK( J ) = WORK( PVT )
                     WORK( PVT ) = DTEMP
                     ITEMP = PIV( PVT )
                     PIV( PVT ) = PIV( J )
                     PIV( J ) = ITEMP
                  END IF
*
                  AJJ = SQRT( AJJ )
                  A( J, J ) = AJJ
*
*                 Compute elements J+1:N of column J.
*
                  IF( J.LT.N ) THEN
                     CALL DGEMV( 'No Trans', N-J, J-K, -ONE, A( J+1, K ), LDA, A( J, K ), LDA, ONE, A( J+1, J ), 1 )
                     CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
                  END IF
*
  170          CONTINUE
*
*              Update trailing matrix, J already incremented
*
               IF( K+JB.LE.N ) THEN
                  CALL DSYRK( 'Lower', 'No Trans', N-J+1, JB, -ONE, A( J, K ), LDA, ONE, A( J, J ), LDA )
               END IF
*
  180       CONTINUE
*
         END IF
      END IF
*
*     Ran to completion, A has full rank
*
      RANK = N
*
      GO TO 200
  190 CONTINUE
*
*     Rank is the number of steps completed.  Set INFO = 1 to signal
*     that the factorization cannot be used to solve a system.
*
      RANK = J - 1
      INFO = 1
*
  200 CONTINUE
      RETURN
*
*     End of DPSTRF
*
      END

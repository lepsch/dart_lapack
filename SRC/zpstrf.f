      SUBROUTINE ZPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   TOL
      int                INFO, LDA, N, RANK
      String             UPLO;
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      DOUBLE PRECISION   WORK( 2*N )
      int                PIV( N )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         ZTEMP
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
      EXTERNAL           ZDSCAL, ZGEMV, ZHERK, ZLACGV, ZPSTF2, ZSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCONJG, MAX, MIN, SQRT, MAXLOC
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
         CALL XERBLA( 'ZPSTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Get block size
*
      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL ZPSTF2( UPLO, N, A( 1, 1 ), LDA, PIV, RANK, TOL, WORK, INFO )
         GO TO 230
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
         DO 110 I = 1, N
            WORK( I ) = DBLE( A( I, I ) )
  110    CONTINUE
         PVT = MAXLOC( WORK( 1:N ), 1 )
         AJJ = DBLE( A( PVT, PVT ) )
         IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
            RANK = 0
            INFO = 1
            GO TO 230
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
*           Compute the Cholesky factorization P**T * A * P = U**H * U
*
            DO 160 K = 1, N, NB
*
*              Account for last block not being NB wide
*
               JB = MIN( NB, N-K+1 )
*
*              Set relevant part of first half of WORK to zero,
*              holds dot products
*
               DO 120 I = K, N
                  WORK( I ) = 0
  120          CONTINUE
*
               DO 150 J = K, K + JB - 1
*
*              Find pivot, test for exit, else swap rows and columns
*              Update dot products, compute possible pivots which are
*              stored in the second half of WORK
*
                  DO 130 I = J, N
*
                     IF( J.GT.K ) THEN
                        WORK( I ) = WORK( I ) + DBLE( DCONJG( A( J-1, I ) )* A( J-1, I ) )
                     END IF
                     WORK( N+I ) = DBLE( A( I, I ) ) - WORK( I )
*
  130             CONTINUE
*
                  IF( J.GT.1 ) THEN
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                     PVT = ITEMP + J - 1
                     AJJ = WORK( N+PVT )
                     IF( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) THEN
                        A( J, J ) = AJJ
                        GO TO 220
                     END IF
                  END IF
*
                  IF( J.NE.PVT ) THEN
*
*                    Pivot OK, so can now swap pivot rows and columns
*
                     A( PVT, PVT ) = A( J, J )
                     CALL ZSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
                     IF( PVT.LT.N ) CALL ZSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA )
                     DO 140 I = J + 1, PVT - 1
                        ZTEMP = DCONJG( A( J, I ) )
                        A( J, I ) = DCONJG( A( I, PVT ) )
                        A( I, PVT ) = ZTEMP
  140                CONTINUE
                     A( J, PVT ) = DCONJG( A( J, PVT ) )
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
                     CALL ZLACGV( J-1, A( 1, J ), 1 )
                     CALL ZGEMV( 'Trans', J-K, N-J, -CONE, A( K, J+1 ), LDA, A( K, J ), 1, CONE, A( J, J+1 ), LDA )
                     CALL ZLACGV( J-1, A( 1, J ), 1 )
                     CALL ZDSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
                  END IF
*
  150          CONTINUE
*
*              Update trailing matrix, J already incremented
*
               IF( K+JB.LE.N ) THEN
                  CALL ZHERK( 'Upper', 'Conj Trans', N-J+1, JB, -ONE, A( K, J ), LDA, ONE, A( J, J ), LDA )
               END IF
*
  160       CONTINUE
*
         ELSE
*
*        Compute the Cholesky factorization P**T * A * P = L * L**H
*
            DO 210 K = 1, N, NB
*
*              Account for last block not being NB wide
*
               JB = MIN( NB, N-K+1 )
*
*              Set relevant part of first half of WORK to zero,
*              holds dot products
*
               DO 170 I = K, N
                  WORK( I ) = 0
  170          CONTINUE
*
               DO 200 J = K, K + JB - 1
*
*              Find pivot, test for exit, else swap rows and columns
*              Update dot products, compute possible pivots which are
*              stored in the second half of WORK
*
                  DO 180 I = J, N
*
                     IF( J.GT.K ) THEN
                        WORK( I ) = WORK( I ) + DBLE( DCONJG( A( I, J-1 ) )* A( I, J-1 ) )
                     END IF
                     WORK( N+I ) = DBLE( A( I, I ) ) - WORK( I )
*
  180             CONTINUE
*
                  IF( J.GT.1 ) THEN
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                     PVT = ITEMP + J - 1
                     AJJ = WORK( N+PVT )
                     IF( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) THEN
                        A( J, J ) = AJJ
                        GO TO 220
                     END IF
                  END IF
*
                  IF( J.NE.PVT ) THEN
*
*                    Pivot OK, so can now swap pivot rows and columns
*
                     A( PVT, PVT ) = A( J, J )
                     CALL ZSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
                     IF( PVT.LT.N ) CALL ZSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
                     DO 190 I = J + 1, PVT - 1
                        ZTEMP = DCONJG( A( I, J ) )
                        A( I, J ) = DCONJG( A( PVT, I ) )
                        A( PVT, I ) = ZTEMP
  190                CONTINUE
                     A( PVT, J ) = DCONJG( A( PVT, J ) )
*
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
                     CALL ZLACGV( J-1, A( J, 1 ), LDA )
                     CALL ZGEMV( 'No Trans', N-J, J-K, -CONE, A( J+1, K ), LDA, A( J, K ), LDA, CONE, A( J+1, J ), 1 )
                     CALL ZLACGV( J-1, A( J, 1 ), LDA )
                     CALL ZDSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
                  END IF
*
  200          CONTINUE
*
*              Update trailing matrix, J already incremented
*
               IF( K+JB.LE.N ) THEN
                  CALL ZHERK( 'Lower', 'No Trans', N-J+1, JB, -ONE, A( J, K ), LDA, ONE, A( J, J ), LDA )
               END IF
*
  210       CONTINUE
*
         END IF
      END IF
*
*     Ran to completion, A has full rank
*
      RANK = N
*
      GO TO 230
  220 CONTINUE
*
*     Rank is the number of steps completed.  Set INFO = 1 to signal
*     that the factorization cannot be used to solve a system.
*
      RANK = J - 1
      INFO = 1
*
  230 CONTINUE
      RETURN
*
*     End of ZPSTRF
*
      END

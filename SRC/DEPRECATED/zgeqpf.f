      SUBROUTINE ZGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, M, N;
*     ..
*     .. Array Arguments ..
      int                JPVT( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I, ITEMP, J, MA, MN, PVT;
      double             TEMP, TEMP2, TOL3Z;
      COMPLEX*16         AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEQR2, ZLARF, ZLARFG, ZSWAP, ZUNM2R
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCMPLX, DCONJG, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DZNRM2;
      EXTERNAL           IDAMAX, DLAMCH, DZNRM2
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQPF', -INFO )
         RETURN
      END IF
*
      MN = MIN( M, N )
      TOL3Z = SQRT(DLAMCH('Epsilon'))
*
*     Move initial columns up front
*
      ITEMP = 1
      DO 10 I = 1, N
         IF( JPVT( I ).NE.0 ) THEN
            IF( I.NE.ITEMP ) THEN
               CALL ZSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            END IF
            ITEMP = ITEMP + 1
         ELSE
            JPVT( I ) = I
         END IF
   10 CONTINUE
      ITEMP = ITEMP - 1
*
*     Compute the QR factorization and update remaining columns
*
      IF( ITEMP.GT.0 ) THEN
         MA = MIN( ITEMP, M )
         CALL ZGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         IF( MA.LT.N ) THEN
            CALL ZUNM2R( 'Left', 'Conjugate transpose', M, N-MA, MA, A, LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO )
         END IF
      END IF
*
      IF( ITEMP.LT.MN ) THEN
*
*        Initialize partial column norms. The first n elements of
*        work store the exact column norms.
*
         DO 20 I = ITEMP + 1, N
            RWORK( I ) = DZNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            RWORK( N+I ) = RWORK( I )
   20    CONTINUE
*
*        Compute factorization
*
         DO 40 I = ITEMP + 1, MN
*
*           Determine ith pivot column and swap if necessary
*
            PVT = ( I-1 ) + IDAMAX( N-I+1, RWORK( I ), 1 )
*
            IF( PVT.NE.I ) THEN
               CALL ZSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               RWORK( PVT ) = RWORK( I )
               RWORK( N+PVT ) = RWORK( N+I )
            END IF
*
*           Generate elementary reflector H(i)
*
            AII = A( I, I )
            CALL ZLARFG( M-I+1, AII, A( MIN( I+1, M ), I ), 1, TAU( I ) )
            A( I, I ) = AII
*
            IF( I.LT.N ) THEN
*
*              Apply H(i) to A(i:m,i+1:n) from the left
*
               AII = A( I, I )
               A( I, I ) = DCMPLX( ONE )
               CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, DCONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK )
               A( I, I ) = AII
            END IF
*
*           Update partial column norms
*
            DO 30 J = I + 1, N
               IF( RWORK( J ).NE.ZERO ) THEN
*
*                 NOTE: The following 4 lines follow from the analysis in
*                 Lapack Working Note 176.
*
                  TEMP = ABS( A( I, J ) ) / RWORK( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( RWORK( J ) / RWORK( N+J ) )**2
                  IF( TEMP2 .LE. TOL3Z ) THEN
                     IF( M-I.GT.0 ) THEN
                        RWORK( J ) = DZNRM2( M-I, A( I+1, J ), 1 )
                        RWORK( N+J ) = RWORK( J )
                     ELSE
                        RWORK( J ) = ZERO
                        RWORK( N+J ) = ZERO
                     END IF
                  ELSE
                     RWORK( J ) = RWORK( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE
*
   40    CONTINUE
      END IF
      RETURN
*
*     End of ZGEQPF
*
      END

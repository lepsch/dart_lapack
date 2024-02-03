      SUBROUTINE DLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2, WI )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB;
      double             SAFMIN, SCALE1, SCALE2, WI, WR1, WR2;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE, TWO;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      double             HALF;
      PARAMETER          ( HALF = ONE / TWO )
      double             FUZZY1;
      PARAMETER          ( FUZZY1 = ONE+1.0D-5 )
*     ..
*     .. Local Scalars ..
      double             A11, A12, A21, A22, ABI22, ANORM, AS11, AS12, AS22, ASCALE, B11, B12, B22, BINV11, BINV22, BMIN, BNORM, BSCALE, BSIZE, C1, C2, C3, C4, C5, DIFF, DISCR, PP, QQ, R, RTMAX, RTMIN, S1, S2, SAFMAX, SHIFT, SS, SUM, WABS, WBIG, WDET, WSCALE, WSIZE, WSMALL;
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
      RTMIN = SQRT( SAFMIN )
      RTMAX = ONE / RTMIN
      SAFMAX = ONE / SAFMIN
*
*     Scale A
*
      ANORM = MAX( ABS( A( 1, 1 ) )+ABS( A( 2, 1 ) ), ABS( A( 1, 2 ) )+ABS( A( 2, 2 ) ), SAFMIN )
      ASCALE = ONE / ANORM
      A11 = ASCALE*A( 1, 1 )
      A21 = ASCALE*A( 2, 1 )
      A12 = ASCALE*A( 1, 2 )
      A22 = ASCALE*A( 2, 2 )
*
*     Perturb B if necessary to insure non-singularity
*
      B11 = B( 1, 1 )
      B12 = B( 1, 2 )
      B22 = B( 2, 2 )
      BMIN = RTMIN*MAX( ABS( B11 ), ABS( B12 ), ABS( B22 ), RTMIN )
      IF( ABS( B11 ).LT.BMIN ) B11 = SIGN( BMIN, B11 )       IF( ABS( B22 ).LT.BMIN ) B22 = SIGN( BMIN, B22 )
*
*     Scale B
*
      BNORM = MAX( ABS( B11 ), ABS( B12 )+ABS( B22 ), SAFMIN )
      BSIZE = MAX( ABS( B11 ), ABS( B22 ) )
      BSCALE = ONE / BSIZE
      B11 = B11*BSCALE
      B12 = B12*BSCALE
      B22 = B22*BSCALE
*
*     Compute larger eigenvalue by method described by C. van Loan
*
*     ( AS is A shifted by -SHIFT*B )
*
      BINV11 = ONE / B11
      BINV22 = ONE / B22
      S1 = A11*BINV11
      S2 = A22*BINV22
      IF( ABS( S1 ).LE.ABS( S2 ) ) THEN
         AS12 = A12 - S1*B12
         AS22 = A22 - S1*B22
         SS = A21*( BINV11*BINV22 )
         ABI22 = AS22*BINV22 - SS*B12
         PP = HALF*ABI22
         SHIFT = S1
      ELSE
         AS12 = A12 - S2*B12
         AS11 = A11 - S2*B11
         SS = A21*( BINV11*BINV22 )
         ABI22 = -SS*B12
         PP = HALF*( AS11*BINV11+ABI22 )
         SHIFT = S2
      END IF
      QQ = SS*AS12
      IF( ABS( PP*RTMIN ).GE.ONE ) THEN
         DISCR = ( RTMIN*PP )**2 + QQ*SAFMIN
         R = SQRT( ABS( DISCR ) )*RTMAX
      ELSE
         IF( PP**2+ABS( QQ ).LE.SAFMIN ) THEN
            DISCR = ( RTMAX*PP )**2 + QQ*SAFMAX
            R = SQRT( ABS( DISCR ) )*RTMIN
         ELSE
            DISCR = PP**2 + QQ
            R = SQRT( ABS( DISCR ) )
         END IF
      END IF
*
*     Note: the test of R in the following IF is to cover the case when
*           DISCR is small and negative and is flushed to zero during
*           the calculation of R.  On machines which have a consistent
*           flush-to-zero threshold and handle numbers above that
*           threshold correctly, it would not be necessary.
*
      IF( DISCR.GE.ZERO .OR. R.EQ.ZERO ) THEN
         SUM = PP + SIGN( R, PP )
         DIFF = PP - SIGN( R, PP )
         WBIG = SHIFT + SUM
*
*        Compute smaller eigenvalue
*
         WSMALL = SHIFT + DIFF
         IF( HALF*ABS( WBIG ).GT.MAX( ABS( WSMALL ), SAFMIN ) ) THEN
            WDET = ( A11*A22-A12*A21 )*( BINV11*BINV22 )
            WSMALL = WDET / WBIG
         END IF
*
*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
*        for WR1.
*
         IF( PP.GT.ABI22 ) THEN
            WR1 = MIN( WBIG, WSMALL )
            WR2 = MAX( WBIG, WSMALL )
         ELSE
            WR1 = MAX( WBIG, WSMALL )
            WR2 = MIN( WBIG, WSMALL )
         END IF
         WI = ZERO
      ELSE
*
*        Complex eigenvalues
*
         WR1 = SHIFT + PP
         WR2 = WR1
         WI = R
      END IF
*
*     Further scaling to avoid underflow and overflow in computing
*     SCALE1 and overflow in computing w*B.
*
*     This scale factor (WSCALE) is bounded from above using C1 and C2,
*     and from below using C3 and C4.
*        C1 implements the condition  s A  must never overflow.
*        C2 implements the condition  w B  must never overflow.
*        C3, with C2,
*           implement the condition that s A - w B must never overflow.
*        C4 implements the condition  s    should not underflow.
*        C5 implements the condition  max(s,|w|) should be at least 2.
*
      C1 = BSIZE*( SAFMIN*MAX( ONE, ASCALE ) )
      C2 = SAFMIN*MAX( ONE, BNORM )
      C3 = BSIZE*SAFMIN
      IF( ASCALE.LE.ONE .AND. BSIZE.LE.ONE ) THEN
         C4 = MIN( ONE, ( ASCALE / SAFMIN )*BSIZE )
      ELSE
         C4 = ONE
      END IF
      IF( ASCALE.LE.ONE .OR. BSIZE.LE.ONE ) THEN
         C5 = MIN( ONE, ASCALE*BSIZE )
      ELSE
         C5 = ONE
      END IF
*
*     Scale first eigenvalue
*
      WABS = ABS( WR1 ) + ABS( WI )
      WSIZE = MAX( SAFMIN, C1, FUZZY1*( WABS*C2+C3 ), MIN( C4, HALF*MAX( WABS, C5 ) ) )
      IF( WSIZE.NE.ONE ) THEN
         WSCALE = ONE / WSIZE
         IF( WSIZE.GT.ONE ) THEN
            SCALE1 = ( MAX( ASCALE, BSIZE )*WSCALE )* MIN( ASCALE, BSIZE )
         ELSE
            SCALE1 = ( MIN( ASCALE, BSIZE )*WSCALE )* MAX( ASCALE, BSIZE )
         END IF
         WR1 = WR1*WSCALE
         IF( WI.NE.ZERO ) THEN
            WI = WI*WSCALE
            WR2 = WR1
            SCALE2 = SCALE1
         END IF
      ELSE
         SCALE1 = ASCALE*BSIZE
         SCALE2 = SCALE1
      END IF
*
*     Scale second eigenvalue (if real)
*
      IF( WI.EQ.ZERO ) THEN
         WSIZE = MAX( SAFMIN, C1, FUZZY1*( ABS( WR2 )*C2+C3 ), MIN( C4, HALF*MAX( ABS( WR2 ), C5 ) ) )
         IF( WSIZE.NE.ONE ) THEN
            WSCALE = ONE / WSIZE
            IF( WSIZE.GT.ONE ) THEN
               SCALE2 = ( MAX( ASCALE, BSIZE )*WSCALE )* MIN( ASCALE, BSIZE )
            ELSE
               SCALE2 = ( MIN( ASCALE, BSIZE )*WSCALE )* MAX( ASCALE, BSIZE )
            END IF
            WR2 = WR2*WSCALE
         ELSE
            SCALE2 = ASCALE*BSIZE
         END IF
      END IF
*
*     End of DLAG2
*
      RETURN
      END

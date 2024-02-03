      void slag2(A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2, WI ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB;
      REAL               SAFMIN, SCALE1, SCALE2, WI, WR1, WR2;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      REAL               HALF;
      const              HALF = ONE / TWO ;
      REAL               FUZZY1;
      const              FUZZY1 = ONE+1.0e-5 ;
      // ..
      // .. Local Scalars ..
      REAL               A11, A12, A21, A22, ABI22, ANORM, AS11, AS12, AS22, ASCALE, B11, B12, B22, BINV11, BINV22, BMIN, BNORM, BSCALE, BSIZE, C1, C2, C3, C4, C5, DIFF, DISCR, PP, QQ, R, RTMAX, RTMIN, S1, S2, SAFMAX, SHIFT, SS, SUM, WABS, WBIG, WDET, WSCALE, WSIZE, WSMALL;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      RTMIN = sqrt( SAFMIN );
      RTMAX = ONE / RTMIN;
      SAFMAX = ONE / SAFMIN;

      // Scale A

      ANORM = max( ( A( 1, 1 ) ).abs()+( A( 2, 1 ) ).abs(), ( A( 1, 2 ) ).abs()+( A( 2, 2 ) ).abs(), SAFMIN );
      ASCALE = ONE / ANORM;
      A11 = ASCALE*A( 1, 1 );
      A21 = ASCALE*A( 2, 1 );
      A12 = ASCALE*A( 1, 2 );
      A22 = ASCALE*A( 2, 2 );

      // Perturb B if necessary to insure non-singularity

      B11 = B( 1, 1 );
      B12 = B( 1, 2 );
      B22 = B( 2, 2 );
      BMIN = RTMIN*max( ( B11 ).abs(), ( B12 ).abs(), ( B22 ).abs(), RTMIN );
      if( ( B11 ).abs() < BMIN ) B11 = SIGN( BMIN, B11 );
      IF( ( B22 ).abs() < BMIN ) B22 = SIGN( BMIN, B22 );

      // Scale B

      BNORM = max( ( B11 ).abs(), ( B12 ).abs()+( B22 ).abs(), SAFMIN );
      BSIZE = max( ( B11 ).abs(), ( B22 ).abs() );
      BSCALE = ONE / BSIZE;
      B11 = B11*BSCALE;
      B12 = B12*BSCALE;
      B22 = B22*BSCALE;

      // Compute larger eigenvalue by method described by C. van Loan

      // ( AS is A shifted by -SHIFT*B )

      BINV11 = ONE / B11;
      BINV22 = ONE / B22;
      S1 = A11*BINV11;
      S2 = A22*BINV22;
      if ( ( S1 ).abs() <= ( S2 ).abs() ) {
         AS12 = A12 - S1*B12;
         AS22 = A22 - S1*B22;
         SS = A21*( BINV11*BINV22 );
         ABI22 = AS22*BINV22 - SS*B12;
         PP = HALF*ABI22;
         SHIFT = S1;
      } else {
         AS12 = A12 - S2*B12;
         AS11 = A11 - S2*B11;
         SS = A21*( BINV11*BINV22 );
         ABI22 = -SS*B12;
         PP = HALF*( AS11*BINV11+ABI22 );
         SHIFT = S2;
      }
      QQ = SS*AS12;
      if ( ( PP*RTMIN ).abs() >= ONE ) {
         DISCR = ( RTMIN*PP )**2 + QQ*SAFMIN;
         R = sqrt( ( DISCR ).abs() )*RTMAX;
      } else {
         if ( PP**2+( QQ ).abs() <= SAFMIN ) {
            DISCR = ( RTMAX*PP )**2 + QQ*SAFMAX;
            R = sqrt( ( DISCR ).abs() )*RTMIN;
         } else {
            DISCR = PP**2 + QQ;
            R = sqrt( ( DISCR ).abs() );
         }
      }

      // Note: the test of R in the following IF is to cover the case when
            // DISCR is small and negative and is flushed to zero during
            // the calculation of R.  On machines which have a consistent
            // flush-to-zero threshold and handle numbers above that
            // threshold correctly, it would not be necessary.

      if ( DISCR >= ZERO || R == ZERO ) {
         SUM = PP + SIGN( R, PP );
         DIFF = PP - SIGN( R, PP );
         WBIG = SHIFT + SUM;

         // Compute smaller eigenvalue

         WSMALL = SHIFT + DIFF;
         if ( HALF*( WBIG ).abs() > max( ( WSMALL ).abs(), SAFMIN ) ) {
            WDET = ( A11*A22-A12*A21 )*( BINV11*BINV22 );
            WSMALL = WDET / WBIG;
         }

         // Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
         // for WR1.

         if ( PP > ABI22 ) {
            WR1 = min( WBIG, WSMALL );
            WR2 = max( WBIG, WSMALL );
         } else {
            WR1 = max( WBIG, WSMALL );
            WR2 = min( WBIG, WSMALL );
         }
         WI = ZERO;
      } else {

         // Complex eigenvalues

         WR1 = SHIFT + PP;
         WR2 = WR1;
         WI = R;
      }

      // Further scaling to avoid underflow and overflow in computing
      // SCALE1 and overflow in computing w*B.

      // This scale factor (WSCALE) is bounded from above using C1 and C2,
      // and from below using C3 and C4.
         // C1 implements the condition  s A  must never overflow.
         // C2 implements the condition  w B  must never overflow.
         // C3, with C2,
            // implement the condition that s A - w B must never overflow.
         // C4 implements the condition  s    should not underflow.
         // C5 implements the condition  max(s,|w|) should be at least 2.

      C1 = BSIZE*( SAFMIN*max( ONE, ASCALE ) );
      C2 = SAFMIN*max( ONE, BNORM );
      C3 = BSIZE*SAFMIN;
      if ( ASCALE <= ONE && BSIZE <= ONE ) {
         C4 = min( ONE, ( ASCALE / SAFMIN )*BSIZE );
      } else {
         C4 = ONE;
      }
      if ( ASCALE <= ONE || BSIZE <= ONE ) {
         C5 = min( ONE, ASCALE*BSIZE );
      } else {
         C5 = ONE;
      }

      // Scale first eigenvalue

      WABS = ( WR1 ).abs() + ( WI ).abs();
      WSIZE = max( SAFMIN, C1, FUZZY1*( WABS*C2+C3 ), min( C4, HALF*max( WABS, C5 ) ) );
      if ( WSIZE != ONE ) {
         WSCALE = ONE / WSIZE;
         if ( WSIZE > ONE ) {
            SCALE1 = ( max( ASCALE, BSIZE )*WSCALE )* min( ASCALE, BSIZE );
         } else {
            SCALE1 = ( min( ASCALE, BSIZE )*WSCALE )* max( ASCALE, BSIZE );
         }
         WR1 = WR1*WSCALE;
         if ( WI != ZERO ) {
            WI = WI*WSCALE;
            WR2 = WR1;
            SCALE2 = SCALE1;
         }
      } else {
         SCALE1 = ASCALE*BSIZE;
         SCALE2 = SCALE1;
      }

      // Scale second eigenvalue (if real)

      if ( WI == ZERO ) {
         WSIZE = max( SAFMIN, C1, FUZZY1*( ( WR2 ).abs()*C2+C3 ), min( C4, HALF*max( ( WR2 ).abs(), C5 ) ) );
         if ( WSIZE != ONE ) {
            WSCALE = ONE / WSIZE;
            if ( WSIZE > ONE ) {
               SCALE2 = ( max( ASCALE, BSIZE )*WSCALE )* min( ASCALE, BSIZE );
            } else {
               SCALE2 = ( min( ASCALE, BSIZE )*WSCALE )* max( ASCALE, BSIZE );
            }
            WR2 = WR2*WSCALE;
         } else {
            SCALE2 = ASCALE*BSIZE;
         }
      }
      return;
      }

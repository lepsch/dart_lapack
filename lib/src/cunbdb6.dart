      void cunbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*);
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ALPHA, REALONE, REALZERO;
      const              ALPHA = 0.83, REALONE = 1.0, REALZERO = 0.0 ;
      COMPLEX            NEGONE, ONE, ZERO;
      const              NEGONE = (-1.0,0.0), ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int                I, IX;
      REAL               EPS, NORM, NORM_NEW, SCL, SSQ;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLASSQ, XERBLA
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0;
      if ( M1 < 0 ) {
         INFO = -1;
      } else if ( M2 < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( INCX1 < 1 ) {
         INFO = -5;
      } else if ( INCX2 < 1 ) {
         INFO = -7;
      } else if ( LDQ1 < max( 1, M1 ) ) {
         INFO = -9;
      } else if ( LDQ2 < max( 1, M2 ) ) {
         INFO = -11;
      } else if ( LWORK < N ) {
         INFO = -13;
      }

      if ( INFO != 0 ) {
         xerbla('CUNBDB6', -INFO );
         return;
      }

      EPS = SLAMCH( 'Precision' );

      // Compute the Euclidean norm of X

      SCL = REALZERO;
      SSQ = REALZERO;
      classq(M1, X1, INCX1, SCL, SSQ );
      classq(M2, X2, INCX2, SCL, SSQ );
      NORM = SCL * sqrt( SSQ );

      // First, project X onto the orthogonal complement of Q's column
      // space

      if ( M1 == 0 ) {
         for (I = 1; I <= N; I++) {
            WORK[I] = ZERO;
         }
      } else {
         cgemv('C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 );
      }

      cgemv('C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 );

      cgemv('N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL CGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 );

      SCL = REALZERO;
      SSQ = REALZERO;
      classq(M1, X1, INCX1, SCL, SSQ );
      classq(M2, X2, INCX2, SCL, SSQ );
      NORM_NEW = SCL * sqrt(SSQ);

      // If projection is sufficiently large in norm, then stop.
      // If projection is zero, then stop.
      // Otherwise, project again.

      if ( NORM_NEW >= ALPHA * NORM ) {
         return;
      }

      if ( NORM_NEW <= N * EPS * NORM ) {
         for (IX = 1; INCX1 < 0 ? IX >= 1 + (M1-1)*INCX1 : IX <= 1 + (M1-1)*INCX1; IX += INCX1) {
           X1[IX] = ZERO;
         }
         for (IX = 1; INCX2 < 0 ? IX >= 1 + (M2-1)*INCX2 : IX <= 1 + (M2-1)*INCX2; IX += INCX2) {
           X2[IX] = ZERO;
         }
         return;
      }

      NORM = NORM_NEW;

      for (I = 1; I <= N; I++) {
         WORK[I] = ZERO;
      }

      if ( M1 == 0 ) {
         for (I = 1; I <= N; I++) {
            WORK[I] = ZERO;
         }
      } else {
         cgemv('C', M1, N, ONE, Q1, LDQ1, X1, INCX1, ZERO, WORK, 1 );
      }

      cgemv('C', M2, N, ONE, Q2, LDQ2, X2, INCX2, ONE, WORK, 1 );

      cgemv('N', M1, N, NEGONE, Q1, LDQ1, WORK, 1, ONE, X1, INCX1 )       CALL CGEMV( 'N', M2, N, NEGONE, Q2, LDQ2, WORK, 1, ONE, X2, INCX2 );

      SCL = REALZERO;
      SSQ = REALZERO;
      classq(M1, X1, INCX1, SCL, SSQ );
      classq(M2, X2, INCX2, SCL, SSQ );
      NORM_NEW = SCL * sqrt(SSQ);

      // If second projection is sufficiently large in norm, then do
      // nothing more. Alternatively, if it shrunk significantly, then
      // truncate it to zero.

      if ( NORM_NEW < ALPHA * NORM ) {
         for (IX = 1; INCX1 < 0 ? IX >= 1 + (M1-1)*INCX1 : IX <= 1 + (M1-1)*INCX1; IX += INCX1) {
            X1[IX] = ZERO;
         }
         for (IX = 1; INCX2 < 0 ? IX >= 1 + (M2-1)*INCX2 : IX <= 1 + (M2-1)*INCX2; IX += INCX2) {
            X2[IX] = ZERO;
         }
      }

      return;
      }

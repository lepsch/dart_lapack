      void slasq2(N, Z, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               Z( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               CBIAS;
      const              CBIAS = 1.50 ;
      REAL               ZERO, HALF, ONE, TWO, FOUR, HUNDRD;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0, TWO = 2.0, FOUR = 4.0, HUNDRD = 100.0 ;
      // ..
      // .. Local Scalars ..
      bool               IEEE;
      int                I0, I4, IINFO, IPN4, ITER, IWHILA, IWHILB, K, KMIN, N0, NBIG, NDIV, NFAIL, PP, SPLT, TTYPE, I1, N1;
      REAL               D, DEE, DEEMIN, DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, E, EMAX, EMIN, EPS, G, OLDEMN, QMAX, QMIN, S, SAFMIN, SIGMA, T, TAU, TEMP, TOL, TOL2, TRACE, ZMAX, TEMPE, TEMPQ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASQ3, SLASRT, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.
      // (in case SLASQ2 is not called by SLASQ1)

      INFO = 0;
      EPS = SLAMCH( 'Precision' );
      SAFMIN = SLAMCH( 'Safe minimum' );
      TOL = EPS*HUNDRD;
      TOL2 = TOL**2;

      if ( N < 0 ) {
         INFO = -1;
         xerbla('SLASQ2', 1 );
         return;
      } else if ( N == 0 ) {
         return;
      } else if ( N == 1 ) {

         // 1-by-1 case.

         if ( Z( 1 ) < ZERO ) {
            INFO = -201;
            xerbla('SLASQ2', 2 );
         }
         return;
      } else if ( N == 2 ) {

         // 2-by-2 case.

         if ( Z( 1 ) < ZERO ) {
            INFO = -201;
            xerbla('SLASQ2', 2 );
            return;
         } else if ( Z( 2 ) < ZERO ) {
            INFO = -202;
            xerbla('SLASQ2', 2 );
            return;
         } else if ( Z( 3 ) < ZERO ) {
           INFO = -203;
           xerbla('SLASQ2', 2 );
           return;
         } else if ( Z( 3 ) > Z( 1 ) ) {
            D = Z( 3 );
            Z[3] = Z( 1 );
            Z[1] = D;
         }
         Z[5] = Z( 1 ) + Z( 2 ) + Z( 3 );
         if ( Z( 2 ) > Z( 3 )*TOL2 ) {
            T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) );
            S = Z( 3 )*( Z( 2 ) / T );
            if ( S <= T ) {
               S = Z( 3 )*( Z( 2 ) / ( T*( ONE+sqrt( ONE+S / T ) ) ) );
            } else {
               S = Z( 3 )*( Z( 2 ) / ( T+sqrt( T )*sqrt( T+S ) ) );
            }
            T = Z( 1 ) + ( S+Z( 2 ) );
            Z[3] = Z( 3 )*( Z( 1 ) / T );
            Z[1] = T;
         }
         Z[2] = Z( 3 );
         Z[6] = Z( 2 ) + Z( 1 );
         return;
      }

      // Check for negative data and compute sums of q's and e's.

      Z[2*N] = ZERO;
      EMIN = Z( 2 );
      QMAX = ZERO;
      ZMAX = ZERO;
      D = ZERO;
      E = ZERO;

      for (K = 1; K <= 2*( N-1 ); K += 2) { // 10
         if ( Z( K ) < ZERO ) {
            INFO = -( 200+K );
            xerbla('SLASQ2', 2 );
            return;
         } else if ( Z( K+1 ) < ZERO ) {
            INFO = -( 200+K+1 );
            xerbla('SLASQ2', 2 );
            return;
         }
         D = D + Z( K );
         E = E + Z( K+1 );
         QMAX = max( QMAX, Z( K ) );
         EMIN = min( EMIN, Z( K+1 ) );
         ZMAX = max( QMAX, ZMAX, Z( K+1 ) );
      } // 10
      if ( Z( 2*N-1 ) < ZERO ) {
         INFO = -( 200+2*N-1 );
         xerbla('SLASQ2', 2 );
         return;
      }
      D = D + Z( 2*N-1 );
      QMAX = max( QMAX, Z( 2*N-1 ) );
      ZMAX = max( QMAX, ZMAX );

      // Check for diagonality.

      if ( E == ZERO ) {
         for (K = 2; K <= N; K++) { // 20
            Z[K] = Z( 2*K-1 );
         } // 20
         slasrt('D', N, Z, IINFO );
         Z[2*N-1] = D;
         return;
      }

      TRACE = D + E;

      // Check for zero data.

      if ( TRACE == ZERO ) {
         Z[2*N-1] = ZERO;
         return;
      }

      // Check whether the machine is IEEE conformable.

      // IEEE = ( ILAENV( 10, 'SLASQ2', 'N', 1, 2, 3, 4 ) == 1 )

      // [11/15/2008] The case IEEE= true has a problem in single precision with
      // some the test matrices of type 16. The double precision code is fine.

      IEEE = false;

      // Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).

      for (K = 2*N; K >= 2; K -= 2) { // 30
         Z[2*K] = ZERO;
         Z[2*K-1] = Z( K );
         Z[2*K-2] = ZERO;
         Z[2*K-3] = Z( K-1 );
      } // 30

      I0 = 1;
      N0 = N;

      // Reverse the qd-array, if warranted.

      if ( CBIAS*Z( 4*I0-3 ) < Z( 4*N0-3 ) ) {
         IPN4 = 4*( I0+N0 );
         for (I4 = 4*I0; I4 <= 2*( I0+N0-1 ); I4 += 4) { // 40
            TEMP = Z( I4-3 );
            Z[I4-3] = Z( IPN4-I4-3 );
            Z[IPN4-I4-3] = TEMP;
            TEMP = Z( I4-1 );
            Z[I4-1] = Z( IPN4-I4-5 );
            Z[IPN4-I4-5] = TEMP;
         } // 40
      }

      // Initial split checking via dqd and Li's test.

      PP = 0;

      for (K = 1; K <= 2; K++) { // 80

         D = Z( 4*N0+PP-3 );
         for (I4 = 4*( N0-1 ) + PP; I4 >= 4*I0 + PP; I4 -= 4) { // 50
            if ( Z( I4-1 ) <= TOL2*D ) {
               Z[I4-1] = -ZERO;
               D = Z( I4-3 );
            } else {
               D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) );
            }
         } // 50

         // dqd maps Z to ZZ plus Li's test.

         EMIN = Z( 4*I0+PP+1 );
         D = Z( 4*I0+PP-3 );
         for (I4 = 4*I0 + PP; I4 <= 4*( N0-1 ) + PP; I4 += 4) { // 60
            Z[I4-2*PP-2] = D + Z( I4-1 );
            if ( Z( I4-1 ) <= TOL2*D ) {
               Z[I4-1] = -ZERO;
               Z[I4-2*PP-2] = D;
               Z[I4-2*PP] = ZERO;
               D = Z( I4+1 );
            } else if ( SAFMIN*Z( I4+1 ) < Z( I4-2*PP-2 ) && SAFMIN*Z( I4-2*PP-2 ) < Z( I4+1 ) ) {
               TEMP = Z( I4+1 ) / Z( I4-2*PP-2 );
               Z[I4-2*PP] = Z( I4-1 )*TEMP;
               D = D*TEMP;
            } else {
               Z[I4-2*PP] = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) );
               D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) );
            }
            EMIN = min( EMIN, Z( I4-2*PP ) );
         } // 60
         Z[4*N0-PP-2] = D;

         // Now find qmax.

         QMAX = Z( 4*I0-PP-2 );
         for (I4 = 4*I0 - PP + 2; 4 < 0 ? I4 >= 4*N0 - PP - 2 : I4 <= 4*N0 - PP - 2; I4 += 4) { // 70
            QMAX = max( QMAX, Z( I4 ) );
         } // 70

         // Prepare for the next iteration on K.

         PP = 1 - PP;
      } // 80

      // Initialise variables to pass to SLASQ3.

      TTYPE = 0;
      DMIN1 = ZERO;
      DMIN2 = ZERO;
      DN    = ZERO;
      DN1   = ZERO;
      DN2   = ZERO;
      G     = ZERO;
      TAU   = ZERO;

      ITER = 2;
      NFAIL = 0;
      NDIV = 2*( N0-I0 );

      for (IWHILA = 1; IWHILA <= N + 1; IWHILA++) { // 160
         if (N0 < 1) GO TO 170;

         // While array unfinished do

         // E(N0) holds the value of SIGMA when submatrix in I0:N0
         // splits from the rest of the array, but is negated.

         DESIG = ZERO;
         if ( N0 == N ) {
            SIGMA = ZERO;
         } else {
            SIGMA = -Z( 4*N0-1 );
         }
         if ( SIGMA < ZERO ) {
            INFO = 1;
            return;
         }

         // Find last unreduced submatrix's top index I0, find QMAX and
         // EMIN. Find Gershgorin-type bound if Q's much greater than E's.

         EMAX = ZERO;
         if ( N0 > I0 ) {
            EMIN = ( Z( 4*N0-5 ) ).abs();
         } else {
            EMIN = ZERO;
         }
         QMIN = Z( 4*N0-3 );
         QMAX = QMIN;
         for (I4 = 4*N0; I4 >= 8; I4 -= 4) { // 90
            if( Z( I4-5 ) <= ZERO ) GO TO 100;
            if ( QMIN >= FOUR*EMAX ) {
               QMIN = min( QMIN, Z( I4-3 ) );
               EMAX = max( EMAX, Z( I4-5 ) );
            }
            QMAX = max( QMAX, Z( I4-7 )+Z( I4-5 ) );
            EMIN = min( EMIN, Z( I4-5 ) );
         } // 90
         I4 = 4;

         } // 100
         I0 = I4 / 4;
         PP = 0;

         if ( N0-I0 > 1 ) {
            DEE = Z( 4*I0-3 );
            DEEMIN = DEE;
            KMIN = I0;
            for (I4 = 4*I0+1; 4 < 0 ? I4 >= 4*N0-3 : I4 <= 4*N0-3; I4 += 4) { // 110
               DEE = Z( I4 )*( DEE /( DEE+Z( I4-2 ) ) );
               if ( DEE <= DEEMIN ) {
                  DEEMIN = DEE;
                  KMIN = ( I4+3 )/4;
               }
            } // 110
            if ( (KMIN-I0)*2 < N0-KMIN && DEEMIN <= HALF*Z(4*N0-3) ) {
               IPN4 = 4*( I0+N0 );
               PP = 2;
               for (I4 = 4*I0; I4 <= 2*( I0+N0-1 ); I4 += 4) { // 120
                  TEMP = Z( I4-3 );
                  Z[I4-3] = Z( IPN4-I4-3 );
                  Z[IPN4-I4-3] = TEMP;
                  TEMP = Z( I4-2 );
                  Z[I4-2] = Z( IPN4-I4-2 );
                  Z[IPN4-I4-2] = TEMP;
                  TEMP = Z( I4-1 );
                  Z[I4-1] = Z( IPN4-I4-5 );
                  Z[IPN4-I4-5] = TEMP;
                  TEMP = Z( I4 );
                  Z[I4] = Z( IPN4-I4-4 );
                  Z[IPN4-I4-4] = TEMP;
               } // 120
            }
         }

         // Put -(initial shift) into DMIN.

         DMIN = -max( ZERO, QMIN-TWO*sqrt( QMIN )*sqrt( EMAX ) );

         // Now I0:N0 is unreduced.
         // PP = 0 for ping, PP = 1 for pong.
         // PP = 2 indicates that flipping was applied to the Z array and
                // and that the tests for deflation upon entry in SLASQ3
                // should not be performed.

         NBIG = 100*( N0-I0+1 );
         for (IWHILB = 1; IWHILB <= NBIG; IWHILB++) { // 140
            if (I0 > N0) GO TO 150;

            // While submatrix unfinished take a good dqds step.

            slasq3(I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, DN2, G, TAU );

            PP = 1 - PP;

            // When EMIN is very small check for splits.

            if ( PP == 0 && N0-I0 >= 3 ) {
               if ( Z( 4*N0 ) <= TOL2*QMAX || Z( 4*N0-1 ) <= TOL2*SIGMA ) {
                  SPLT = I0 - 1;
                  QMAX = Z( 4*I0-3 );
                  EMIN = Z( 4*I0-1 );
                  OLDEMN = Z( 4*I0 );
                  for (I4 = 4*I0; I4 <= 4*( N0-3 ); I4 += 4) { // 130
                     if ( Z( I4 ) <= TOL2*Z( I4-3 ) || Z( I4-1 ) <= TOL2*SIGMA ) {
                        Z[I4-1] = -SIGMA;
                        SPLT = I4 / 4;
                        QMAX = ZERO;
                        EMIN = Z( I4+3 );
                        OLDEMN = Z( I4+4 );
                     } else {
                        QMAX = max( QMAX, Z( I4+1 ) );
                        EMIN = min( EMIN, Z( I4-1 ) );
                        OLDEMN = min( OLDEMN, Z( I4 ) );
                     }
                  } // 130
                  Z[4*N0-1] = EMIN;
                  Z[4*N0] = OLDEMN;
                  I0 = SPLT + 1;
               }
            }

         } // 140

         INFO = 2;

         // Maximum number of iterations exceeded, restore the shift
         // SIGMA and place the new d's and e's in a qd array.
         // This might need to be done for several blocks

         I1 = I0;
         N1 = N0;
         } // 145
         TEMPQ = Z( 4*I0-3 );
         Z[4*I0-3] = Z( 4*I0-3 ) + SIGMA;
         for (K = I0+1; K <= N0; K++) {
            TEMPE = Z( 4*K-5 );
            Z[4*K-5] = Z( 4*K-5 ) * (TEMPQ / Z( 4*K-7 ));
            TEMPQ = Z( 4*K-3 );
            Z[4*K-3] = Z( 4*K-3 ) + SIGMA + TEMPE - Z( 4*K-5 );
         }

         // Prepare to do this on the previous block if there is one

         if ( I1 > 1 ) {
            N1 = I1-1;
            while (( I1 >= 2 ) && ( Z(4*I1-5) >= ZERO )) {
               I1 = I1 - 1;
            }
            if ( I1 >= 1 ) {
               SIGMA = -Z(4*N1-1);
               GO TO 145;
            }
         }

         for (K = 1; K <= N; K++) {
            Z[2*K-1] = Z( 4*K-3 );

         // Only the block 1..N0 is unfinished.  The rest of the e's
         // must be essentially zero, although sometimes other data
         // has been stored in them.

            if ( K < N0 ) {
               Z[2*K] = Z( 4*K-1 );
            } else {
               Z[2*K] = 0;
            }
         }
         return;

         // end IWHILB

         } // 150

      } // 160

      INFO = 3;
      return;

      // end IWHILA

      } // 170

      // Move q's to the front.

      for (K = 2; K <= N; K++) { // 180
         Z[K] = Z( 4*K-3 );
      } // 180

      // Sort and compute sum of eigenvalues.

      slasrt('D', N, Z, IINFO );

      E = ZERO;
      for (K = N; K >= 1; K--) { // 190
         E = E + Z( K );
      } // 190

      // Store trace, sum(eigenvalues) and information on performance.

      Z[2*N+1] = TRACE;
      Z[2*N+2] = E;
      Z[2*N+3] = REAL( ITER );
      Z[2*N+4] = REAL( NDIV ) / REAL( N**2 );
      Z[2*N+5] = HUNDRD*NFAIL / REAL( ITER );
      return;
      }

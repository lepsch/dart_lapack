      SUBROUTINE SLASQ2( N, Z, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               Z( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               CBIAS
      const              CBIAS = 1.50E0 ;
      REAL               ZERO, HALF, ONE, TWO, FOUR, HUNDRD
      const              ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0, TWO = 2.0E0, FOUR = 4.0E0, HUNDRD = 100.0E0 ;
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
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.
      // (in case SLASQ2 is not called by SLASQ1)

      INFO = 0
      EPS = SLAMCH( 'Precision' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2

      if ( N.LT.0 ) {
         INFO = -1
         CALL XERBLA( 'SLASQ2', 1 )
         RETURN
      } else if ( N.EQ.0 ) {
         RETURN
      } else if ( N.EQ.1 ) {

         // 1-by-1 case.

         if ( Z( 1 ).LT.ZERO ) {
            INFO = -201
            CALL XERBLA( 'SLASQ2', 2 )
         }
         RETURN
      } else if ( N.EQ.2 ) {

         // 2-by-2 case.

         if ( Z( 1 ).LT.ZERO ) {
            INFO = -201
            CALL XERBLA( 'SLASQ2', 2 )
            RETURN
         } else if ( Z( 2 ).LT.ZERO ) {
            INFO = -202
            CALL XERBLA( 'SLASQ2', 2 )
            RETURN
         } else if ( Z( 3 ).LT.ZERO ) {
           INFO = -203
           CALL XERBLA( 'SLASQ2', 2 )
           RETURN
         } else if ( Z( 3 ).GT.Z( 1 ) ) {
            D = Z( 3 )
            Z( 3 ) = Z( 1 )
            Z( 1 ) = D
         }
         Z( 5 ) = Z( 1 ) + Z( 2 ) + Z( 3 )
         if ( Z( 2 ).GT.Z( 3 )*TOL2 ) {
            T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) )
            S = Z( 3 )*( Z( 2 ) / T )
            if ( S.LE.T ) {
               S = Z( 3 )*( Z( 2 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
            } else {
               S = Z( 3 )*( Z( 2 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
            }
            T = Z( 1 ) + ( S+Z( 2 ) )
            Z( 3 ) = Z( 3 )*( Z( 1 ) / T )
            Z( 1 ) = T
         }
         Z( 2 ) = Z( 3 )
         Z( 6 ) = Z( 2 ) + Z( 1 )
         RETURN
      }

      // Check for negative data and compute sums of q's and e's.

      Z( 2*N ) = ZERO
      EMIN = Z( 2 )
      QMAX = ZERO
      ZMAX = ZERO
      D = ZERO
      E = ZERO

      DO 10 K = 1, 2*( N-1 ), 2
         if ( Z( K ).LT.ZERO ) {
            INFO = -( 200+K )
            CALL XERBLA( 'SLASQ2', 2 )
            RETURN
         } else if ( Z( K+1 ).LT.ZERO ) {
            INFO = -( 200+K+1 )
            CALL XERBLA( 'SLASQ2', 2 )
            RETURN
         }
         D = D + Z( K )
         E = E + Z( K+1 )
         QMAX = MAX( QMAX, Z( K ) )
         EMIN = MIN( EMIN, Z( K+1 ) )
         ZMAX = MAX( QMAX, ZMAX, Z( K+1 ) )
   10 CONTINUE
      if ( Z( 2*N-1 ).LT.ZERO ) {
         INFO = -( 200+2*N-1 )
         CALL XERBLA( 'SLASQ2', 2 )
         RETURN
      }
      D = D + Z( 2*N-1 )
      QMAX = MAX( QMAX, Z( 2*N-1 ) )
      ZMAX = MAX( QMAX, ZMAX )

      // Check for diagonality.

      if ( E.EQ.ZERO ) {
         DO 20 K = 2, N
            Z( K ) = Z( 2*K-1 )
   20    CONTINUE
         CALL SLASRT( 'D', N, Z, IINFO )
         Z( 2*N-1 ) = D
         RETURN
      }

      TRACE = D + E

      // Check for zero data.

      if ( TRACE.EQ.ZERO ) {
         Z( 2*N-1 ) = ZERO
         RETURN
      }

      // Check whether the machine is IEEE conformable.

      // IEEE = ( ILAENV( 10, 'SLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 )

      // [11/15/2008] The case IEEE=.TRUE. has a problem in single precision with
      // some the test matrices of type 16. The double precision code is fine.

      IEEE = .FALSE.

      // Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).

      DO 30 K = 2*N, 2, -2
         Z( 2*K ) = ZERO
         Z( 2*K-1 ) = Z( K )
         Z( 2*K-2 ) = ZERO
         Z( 2*K-3 ) = Z( K-1 )
   30 CONTINUE

      I0 = 1
      N0 = N

      // Reverse the qd-array, if warranted.

      if ( CBIAS*Z( 4*I0-3 ).LT.Z( 4*N0-3 ) ) {
         IPN4 = 4*( I0+N0 )
         DO 40 I4 = 4*I0, 2*( I0+N0-1 ), 4
            TEMP = Z( I4-3 )
            Z( I4-3 ) = Z( IPN4-I4-3 )
            Z( IPN4-I4-3 ) = TEMP
            TEMP = Z( I4-1 )
            Z( I4-1 ) = Z( IPN4-I4-5 )
            Z( IPN4-I4-5 ) = TEMP
   40    CONTINUE
      }

      // Initial split checking via dqd and Li's test.

      PP = 0

      DO 80 K = 1, 2

         D = Z( 4*N0+PP-3 )
         DO 50 I4 = 4*( N0-1 ) + PP, 4*I0 + PP, -4
            if ( Z( I4-1 ).LE.TOL2*D ) {
               Z( I4-1 ) = -ZERO
               D = Z( I4-3 )
            } else {
               D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) )
            }
   50    CONTINUE

         // dqd maps Z to ZZ plus Li's test.

         EMIN = Z( 4*I0+PP+1 )
         D = Z( 4*I0+PP-3 )
         DO 60 I4 = 4*I0 + PP, 4*( N0-1 ) + PP, 4
            Z( I4-2*PP-2 ) = D + Z( I4-1 )
            if ( Z( I4-1 ).LE.TOL2*D ) {
               Z( I4-1 ) = -ZERO
               Z( I4-2*PP-2 ) = D
               Z( I4-2*PP ) = ZERO
               D = Z( I4+1 )
            } else if ( SAFMIN*Z( I4+1 ).LT.Z( I4-2*PP-2 ) .AND. SAFMIN*Z( I4-2*PP-2 ).LT.Z( I4+1 ) ) {
               TEMP = Z( I4+1 ) / Z( I4-2*PP-2 )
               Z( I4-2*PP ) = Z( I4-1 )*TEMP
               D = D*TEMP
            } else {
               Z( I4-2*PP ) = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) )
               D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) )
            }
            EMIN = MIN( EMIN, Z( I4-2*PP ) )
   60    CONTINUE
         Z( 4*N0-PP-2 ) = D

         // Now find qmax.

         QMAX = Z( 4*I0-PP-2 )
         DO 70 I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
            QMAX = MAX( QMAX, Z( I4 ) )
   70    CONTINUE

         // Prepare for the next iteration on K.

         PP = 1 - PP
   80 CONTINUE

      // Initialise variables to pass to SLASQ3.

      TTYPE = 0
      DMIN1 = ZERO
      DMIN2 = ZERO
      DN    = ZERO
      DN1   = ZERO
      DN2   = ZERO
      G     = ZERO
      TAU   = ZERO

      ITER = 2
      NFAIL = 0
      NDIV = 2*( N0-I0 )

      DO 160 IWHILA = 1, N + 1
         IF( N0.LT.1 ) GO TO 170

         // While array unfinished do

         // E(N0) holds the value of SIGMA when submatrix in I0:N0
         // splits from the rest of the array, but is negated.

         DESIG = ZERO
         if ( N0.EQ.N ) {
            SIGMA = ZERO
         } else {
            SIGMA = -Z( 4*N0-1 )
         }
         if ( SIGMA.LT.ZERO ) {
            INFO = 1
            RETURN
         }

         // Find last unreduced submatrix's top index I0, find QMAX and
         // EMIN. Find Gershgorin-type bound if Q's much greater than E's.

         EMAX = ZERO
         if ( N0.GT.I0 ) {
            EMIN = ABS( Z( 4*N0-5 ) )
         } else {
            EMIN = ZERO
         }
         QMIN = Z( 4*N0-3 )
         QMAX = QMIN
         DO 90 I4 = 4*N0, 8, -4
            IF( Z( I4-5 ).LE.ZERO ) GO TO 100
            if ( QMIN.GE.FOUR*EMAX ) {
               QMIN = MIN( QMIN, Z( I4-3 ) )
               EMAX = MAX( EMAX, Z( I4-5 ) )
            }
            QMAX = MAX( QMAX, Z( I4-7 )+Z( I4-5 ) )
            EMIN = MIN( EMIN, Z( I4-5 ) )
   90    CONTINUE
         I4 = 4

  100    CONTINUE
         I0 = I4 / 4
         PP = 0

         if ( N0-I0.GT.1 ) {
            DEE = Z( 4*I0-3 )
            DEEMIN = DEE
            KMIN = I0
            DO 110 I4 = 4*I0+1, 4*N0-3, 4
               DEE = Z( I4 )*( DEE /( DEE+Z( I4-2 ) ) )
               if ( DEE.LE.DEEMIN ) {
                  DEEMIN = DEE
                  KMIN = ( I4+3 )/4
               }
  110       CONTINUE
            if ( (KMIN-I0)*2.LT.N0-KMIN .AND. DEEMIN.LE.HALF*Z(4*N0-3) ) {
               IPN4 = 4*( I0+N0 )
               PP = 2
               DO 120 I4 = 4*I0, 2*( I0+N0-1 ), 4
                  TEMP = Z( I4-3 )
                  Z( I4-3 ) = Z( IPN4-I4-3 )
                  Z( IPN4-I4-3 ) = TEMP
                  TEMP = Z( I4-2 )
                  Z( I4-2 ) = Z( IPN4-I4-2 )
                  Z( IPN4-I4-2 ) = TEMP
                  TEMP = Z( I4-1 )
                  Z( I4-1 ) = Z( IPN4-I4-5 )
                  Z( IPN4-I4-5 ) = TEMP
                  TEMP = Z( I4 )
                  Z( I4 ) = Z( IPN4-I4-4 )
                  Z( IPN4-I4-4 ) = TEMP
  120          CONTINUE
            }
         }

         // Put -(initial shift) into DMIN.

         DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN )*SQRT( EMAX ) )

         // Now I0:N0 is unreduced.
         // PP = 0 for ping, PP = 1 for pong.
         // PP = 2 indicates that flipping was applied to the Z array and
                // and that the tests for deflation upon entry in SLASQ3
                // should not be performed.

         NBIG = 100*( N0-I0+1 )
         DO 140 IWHILB = 1, NBIG
            IF( I0.GT.N0 ) GO TO 150

            // While submatrix unfinished take a good dqds step.

            CALL SLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, DN2, G, TAU )

            PP = 1 - PP

            // When EMIN is very small check for splits.

            if ( PP.EQ.0 .AND. N0-I0.GE.3 ) {
               if ( Z( 4*N0 ).LE.TOL2*QMAX .OR. Z( 4*N0-1 ).LE.TOL2*SIGMA ) {
                  SPLT = I0 - 1
                  QMAX = Z( 4*I0-3 )
                  EMIN = Z( 4*I0-1 )
                  OLDEMN = Z( 4*I0 )
                  DO 130 I4 = 4*I0, 4*( N0-3 ), 4
                     if ( Z( I4 ).LE.TOL2*Z( I4-3 ) .OR. Z( I4-1 ).LE.TOL2*SIGMA ) {
                        Z( I4-1 ) = -SIGMA
                        SPLT = I4 / 4
                        QMAX = ZERO
                        EMIN = Z( I4+3 )
                        OLDEMN = Z( I4+4 )
                     } else {
                        QMAX = MAX( QMAX, Z( I4+1 ) )
                        EMIN = MIN( EMIN, Z( I4-1 ) )
                        OLDEMN = MIN( OLDEMN, Z( I4 ) )
                     }
  130             CONTINUE
                  Z( 4*N0-1 ) = EMIN
                  Z( 4*N0 ) = OLDEMN
                  I0 = SPLT + 1
               }
            }

  140    CONTINUE

         INFO = 2

         // Maximum number of iterations exceeded, restore the shift
         // SIGMA and place the new d's and e's in a qd array.
         // This might need to be done for several blocks

         I1 = I0
         N1 = N0
 145     CONTINUE
         TEMPQ = Z( 4*I0-3 )
         Z( 4*I0-3 ) = Z( 4*I0-3 ) + SIGMA
         DO K = I0+1, N0
            TEMPE = Z( 4*K-5 )
            Z( 4*K-5 ) = Z( 4*K-5 ) * (TEMPQ / Z( 4*K-7 ))
            TEMPQ = Z( 4*K-3 )
            Z( 4*K-3 ) = Z( 4*K-3 ) + SIGMA + TEMPE - Z( 4*K-5 )
         END DO

         // Prepare to do this on the previous block if there is one

         if ( I1.GT.1 ) {
            N1 = I1-1
            DO WHILE( ( I1.GE.2 ) .AND. ( Z(4*I1-5).GE.ZERO ) )
               I1 = I1 - 1
            END DO
            if ( I1.GE.1 ) {
               SIGMA = -Z(4*N1-1)
               GO TO 145
            }
         }

         DO K = 1, N
            Z( 2*K-1 ) = Z( 4*K-3 )

         // Only the block 1..N0 is unfinished.  The rest of the e's
         // must be essentially zero, although sometimes other data
         // has been stored in them.

            if ( K.LT.N0 ) {
               Z( 2*K ) = Z( 4*K-1 )
            } else {
               Z( 2*K ) = 0
            }
         END DO
         RETURN

         // end IWHILB

  150    CONTINUE

  160 CONTINUE

      INFO = 3
      RETURN

      // end IWHILA

  170 CONTINUE

      // Move q's to the front.

      DO 180 K = 2, N
         Z( K ) = Z( 4*K-3 )
  180 CONTINUE

      // Sort and compute sum of eigenvalues.

      CALL SLASRT( 'D', N, Z, IINFO )

      E = ZERO
      DO 190 K = N, 1, -1
         E = E + Z( K )
  190 CONTINUE

      // Store trace, sum(eigenvalues) and information on performance.

      Z( 2*N+1 ) = TRACE
      Z( 2*N+2 ) = E
      Z( 2*N+3 ) = REAL( ITER )
      Z( 2*N+4 ) = REAL( NDIV ) / REAL( N**2 )
      Z( 2*N+5 ) = HUNDRD*NFAIL / REAL( ITER )
      RETURN

      // End of SLASQ2

      }

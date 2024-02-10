      void shgeqz(final int JOB, final int COMPQ, final int COMPZ, final int N, final int ILO, final int IHI, final Matrix<double> H, final int LDH, final Matrix<double> T, final int LDT, final int ALPHAR, final int ALPHAI, final int BETA, final Matrix<double> Q, final int LDQ, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             COMPQ, COMPZ, JOB;
      int                IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N;
      double               ALPHAI( * ), ALPHAR( * ), BETA( * ), H( LDH, * ), Q( LDQ, * ), T( LDT, * ), WORK( * ), Z( LDZ, * );
      // ..

// $                     SAFETY = 1.0 )
      double               HALF, ZERO, ONE, SAFETY;
      const              HALF = 0.5, ZERO = 0.0, ONE = 1.0, SAFETY = 1.0e+2 ;
      bool               ILAZR2, ILAZRO, ILPIVT, ILQ, ILSCHR, ILZ, LQUERY;
      int                ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST, ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER, JR, MAXIT;
      double               A11, A12, A1I, A1R, A21, A22, A2I, A2R, AD11, AD11L, AD12, AD12L, AD21, AD21L, AD22, AD22L, AD32L, AN, ANORM, ASCALE, ATOL, B11, B1A, B1I, B1R, B22, B2A, B2I, B2R, BN, BNORM, BSCALE, BTOL, C, C11I, C11R, C12, C21, C22I, C22R, CL, CQ, CR, CZ, ESHIFT, S, S1, S1INV, S2, SAFMAX, SAFMIN, SCALE, SL, SQI, SQR, SR, SZI, SZR, T1, T2, T3, TAU, TEMP, TEMP2, TEMPI, TEMPR, U1, U12, U12L, U2, ULP, VS, W11, W12, W21, W22, WABS, WI, WR, WR2;
      double               V( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANHS, SLAPY2, SLAPY3, SROUNDUP_LWORK;
      // EXTERNAL lsame, SLAMCH, SLANHS, SLAPY2, SLAPY3, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAG2, SLARFG, SLARTG, SLASET, SLASV2, SROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT

      // Decode JOB, COMPQ, COMPZ

      if ( lsame( JOB, 'E' ) ) {
         ILSCHR = false;
         ISCHUR = 1;
      } else if ( lsame( JOB, 'S' ) ) {
         ILSCHR = true;
         ISCHUR = 2;
      } else {
         ISCHUR = 0;
      }

      if ( lsame( COMPQ, 'N' ) ) {
         ILQ = false;
         ICOMPQ = 1;
      } else if ( lsame( COMPQ, 'V' ) ) {
         ILQ = true;
         ICOMPQ = 2;
      } else if ( lsame( COMPQ, 'I' ) ) {
         ILQ = true;
         ICOMPQ = 3;
      } else {
         ICOMPQ = 0;
      }

      if ( lsame( COMPZ, 'N' ) ) {
         ILZ = false;
         ICOMPZ = 1;
      } else if ( lsame( COMPZ, 'V' ) ) {
         ILZ = true;
         ICOMPZ = 2;
      } else if ( lsame( COMPZ, 'I' ) ) {
         ILZ = true;
         ICOMPZ = 3;
      } else {
         ICOMPZ = 0;
      }

      // Check Argument Values

      INFO = 0;
      WORK[1] = max( 1, N );
      LQUERY = ( LWORK == -1 );
      if ( ISCHUR == 0 ) {
         INFO = -1;
      } else if ( ICOMPQ == 0 ) {
         INFO = -2;
      } else if ( ICOMPZ == 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( ILO < 1 ) {
         INFO = -5;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO = -6;
      } else if ( LDH < N ) {
         INFO = -8;
      } else if ( LDT < N ) {
         INFO = -10;
      } else if ( LDQ < 1 || ( ILQ && LDQ < N ) ) {
         INFO = -15;
      } else if ( LDZ < 1 || ( ILZ && LDZ < N ) ) {
         INFO = -17;
      } else if ( LWORK < max( 1, N ) && !LQUERY ) {
         INFO = -19;
      }
      if ( INFO != 0 ) {
         xerbla('SHGEQZ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N <= 0 ) {
         WORK[1] = double( 1 );
         return;
      }

      // Initialize Q and Z

      if (ICOMPQ == 3) slaset( 'Full', N, N, ZERO, ONE, Q, LDQ );
      IF( ICOMPZ == 3 ) slaset( 'Full', N, N, ZERO, ONE, Z, LDZ );

      // Machine Constants

      IN = IHI + 1 - ILO;
      SAFMIN = SLAMCH( 'S' );
      SAFMAX = ONE / SAFMIN;
      ULP = SLAMCH( 'E' )*SLAMCH( 'B' );
      ANORM = SLANHS( 'F', IN, H( ILO, ILO ), LDH, WORK );
      BNORM = SLANHS( 'F', IN, T( ILO, ILO ), LDT, WORK );
      ATOL = max( SAFMIN, ULP*ANORM );
      BTOL = max( SAFMIN, ULP*BNORM );
      ASCALE = ONE / max( SAFMIN, ANORM );
      BSCALE = ONE / max( SAFMIN, BNORM );

      // Set Eigenvalues IHI+1:N

      for (J = IHI + 1; J <= N; J++) { // 30
         if ( T( J, J ) < ZERO ) {
            if ( ILSCHR ) {
               for (JR = 1; JR <= J; JR++) { // 10
                  H[JR][J] = -H( JR, J );
                  T[JR][J] = -T( JR, J );
               } // 10
            } else {
               H[J][J] = -H( J, J );
               T[J][J] = -T( J, J );
            }
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 20
                  Z[JR][J] = -Z( JR, J );
               } // 20
            }
         }
         ALPHAR[J] = H( J, J );
         ALPHAI[J] = ZERO;
         BETA[J] = T( J, J );
      } // 30

      // If IHI < ILO, skip QZ steps

      if (IHI < ILO) GO TO 380;

      // MAIN QZ ITERATION LOOP

      // Initialize dynamic indices

      // Eigenvalues ILAST+1:N have been found.
      //    Column operations modify rows IFRSTM:whatever.
      //    Row operations modify columns whatever:ILASTM.

      // If only eigenvalues are being computed, then
      //    IFRSTM is the row of the last splitting row above row ILAST;
      //    this is always at least ILO.
      // IITER counts iterations since the last eigenvalue was found,
      //    to tell when to use an extraordinary shift.
      // MAXIT is the maximum number of QZ sweeps allowed.

      ILAST = IHI;
      if ( ILSCHR ) {
         IFRSTM = 1;
         ILASTM = N;
      } else {
         IFRSTM = ILO;
         ILASTM = IHI;
      }
      IITER = 0;
      ESHIFT = ZERO;
      MAXIT = 30*( IHI-ILO+1 );

      for (JITER = 1; JITER <= MAXIT; JITER++) { // 360

         // Split the matrix if possible.

         // Two tests:
         //    1: H(j,j-1)=0  or  j=ILO
         //    2: T(j,j)=0

         if ( ILAST == ILO ) {

            // Special case: j=ILAST

            GO TO 80;
         } else {
            if ( ( H( ILAST, ILAST-1 ) ).abs() <= max( SAFMIN, ULP*(  ( H( ILAST, ILAST ) ).abs() + ( H( ILAST-1, ILAST-1 ) ).abs() ) ) ) {
               H[ILAST][ILAST-1] = ZERO;
               GO TO 80;
            }
         }

         if ( ( T( ILAST, ILAST ) ).abs() <= BTOL ) {
            T[ILAST][ILAST] = ZERO;
            GO TO 70;
         }

         // General case: j<ILAST

         for (J = ILAST - 1; J >= ILO; J--) { // 60

            // Test 1: for H(j,j-1)=0 or j=ILO

            if ( J == ILO ) {
               ILAZRO = true;
            } else {
               if ( ( H( J, J-1 ) ).abs() <= max( SAFMIN, ULP*(  ( H( J, J ) ).abs() + ( H( J-1, J-1 ) ).abs() ) ) ) {
                  H[J][J-1] = ZERO;
                  ILAZRO = true;
               } else {
                  ILAZRO = false;
               }
            }

            // Test 2: for T(j,j)=0

            if ( ( T( J, J ) ).abs() < BTOL ) {
               T[J][J] = ZERO;

               // Test 1a: Check for 2 consecutive small subdiagonals in A

               ILAZR2 = false;
               if ( !ILAZRO ) {
                  TEMP = ( H( J, J-1 ) ).abs();
                  TEMP2 = ( H( J, J ) ).abs();
                  TEMPR = max( TEMP, TEMP2 );
                  if ( TEMPR < ONE && TEMPR != ZERO ) {
                     TEMP = TEMP / TEMPR;
                     TEMP2 = TEMP2 / TEMPR;
                  }
                  if( TEMP*( ASCALE*( H( J+1, J ) ).abs() ) <= TEMP2* ( ASCALE*ATOL ) )ILAZR2 = true;
               }

               // If both tests pass (1 & 2), i.e., the leading diagonal
               // element of B in the block is zero, split a 1x1 block off
               // at the top. (I.e., at the J-th row/column) The leading
               // diagonal element of the remainder can also be zero, so
               // this may have to be done repeatedly.

               if ( ILAZRO || ILAZR2 ) {
                  for (JCH = J; JCH <= ILAST - 1; JCH++) { // 40
                     TEMP = H( JCH, JCH );
                     slartg(TEMP, H( JCH+1, JCH ), C, S, H( JCH, JCH ) );
                     H[JCH+1][JCH] = ZERO;
                     srot(ILASTM-JCH, H( JCH, JCH+1 ), LDH, H( JCH+1, JCH+1 ), LDH, C, S );
                     srot(ILASTM-JCH, T( JCH, JCH+1 ), LDT, T( JCH+1, JCH+1 ), LDT, C, S )                      IF( ILQ ) CALL SROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1, C, S );
                     if (ILAZR2) H( JCH, JCH-1 ) = H( JCH, JCH-1 )*C;
                     ILAZR2 = false;
                     if ( ( T( JCH+1, JCH+1 ) ).abs() >= BTOL ) {
                        if ( JCH+1 >= ILAST ) {
                           GO TO 80;
                        } else {
                           IFIRST = JCH + 1;
                           GO TO 110;
                        }
                     }
                     T[JCH+1][JCH+1] = ZERO;
                  } // 40
                  GO TO 70;
               } else {

                  // Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                  // Then process as in the case T(ILAST,ILAST)=0

                  for (JCH = J; JCH <= ILAST - 1; JCH++) { // 50
                     TEMP = T( JCH, JCH+1 );
                     slartg(TEMP, T( JCH+1, JCH+1 ), C, S, T( JCH, JCH+1 ) );
                     T[JCH+1][JCH+1] = ZERO;
                     if (JCH < ILASTM-1) srot( ILASTM-JCH-1, T( JCH, JCH+2 ), LDT, T( JCH+1, JCH+2 ), LDT, C, S );
                     srot(ILASTM-JCH+2, H( JCH, JCH-1 ), LDH, H( JCH+1, JCH-1 ), LDH, C, S )                      IF( ILQ ) CALL SROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1, C, S );
                     TEMP = H( JCH+1, JCH );
                     slartg(TEMP, H( JCH+1, JCH-1 ), C, S, H( JCH+1, JCH ) );
                     H[JCH+1][JCH-1] = ZERO;
                     srot(JCH+1-IFRSTM, H( IFRSTM, JCH ), 1, H( IFRSTM, JCH-1 ), 1, C, S );
                     srot(JCH-IFRSTM, T( IFRSTM, JCH ), 1, T( IFRSTM, JCH-1 ), 1, C, S )                      IF( ILZ ) CALL SROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1, C, S );
                  } // 50
                  GO TO 70;
               }
            } else if ( ILAZRO ) {

               // Only test 1 passed -- work on J:ILAST

               IFIRST = J;
               GO TO 110;
            }

            // Neither test passed -- try next J

         } // 60

         // (Drop-through is "impossible")

         INFO = N + 1;
         GO TO 420;

         // T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
         // 1x1 block.

         } // 70
         TEMP = H( ILAST, ILAST );
         slartg(TEMP, H( ILAST, ILAST-1 ), C, S, H( ILAST, ILAST ) );
         H[ILAST][ILAST-1] = ZERO;
         srot(ILAST-IFRSTM, H( IFRSTM, ILAST ), 1, H( IFRSTM, ILAST-1 ), 1, C, S );
         srot(ILAST-IFRSTM, T( IFRSTM, ILAST ), 1, T( IFRSTM, ILAST-1 ), 1, C, S )          IF( ILZ ) CALL SROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S );

         // H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
         //                       and BETA

         } // 80
         if ( T( ILAST, ILAST ) < ZERO ) {
            if ( ILSCHR ) {
               for (J = IFRSTM; J <= ILAST; J++) { // 90
                  H[J][ILAST] = -H( J, ILAST );
                  T[J][ILAST] = -T( J, ILAST );
               } // 90
            } else {
               H[ILAST][ILAST] = -H( ILAST, ILAST );
               T[ILAST][ILAST] = -T( ILAST, ILAST );
            }
            if ( ILZ ) {
               for (J = 1; J <= N; J++) { // 100
                  Z[J][ILAST] = -Z( J, ILAST );
               } // 100
            }
         }
         ALPHAR[ILAST] = H( ILAST, ILAST );
         ALPHAI[ILAST] = ZERO;
         BETA[ILAST] = T( ILAST, ILAST );

         // Go to next block -- exit if finished.

         ILAST = ILAST - 1;
         if (ILAST < ILO) GO TO 380;

         // Reset counters

         IITER = 0;
         ESHIFT = ZERO;
         if ( !ILSCHR ) {
            ILASTM = ILAST;
            if (IFRSTM > ILAST) IFRSTM = ILO;
         }
         GO TO 350;

         // QZ step

         // This iteration only involves rows/columns IFIRST:ILAST. We
         // assume IFIRST < ILAST, and that the diagonal of B is non-zero.

         } // 110
         IITER = IITER + 1;
         if ( !ILSCHR ) {
            IFRSTM = IFIRST;
         }

         // Compute single shifts.

         // At this point, IFIRST < ILAST, and the diagonal elements of
         // T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
         // magnitude)

         if ( ( IITER ~/ 10 )*10 == IITER ) {

            // Exceptional shift.  Chosen for no particularly good reason.
            // (Single shift only.)

            if( ( double( MAXIT )*SAFMIN )*( H( ILAST, ILAST-1 ) ).abs() < ( T( ILAST-1, ILAST-1 ) ).abs() ) {
               ESHIFT = H( ILAST, ILAST-1 ) / T( ILAST-1, ILAST-1 );
            } else {
               ESHIFT = ESHIFT + ONE / ( SAFMIN*REAL( MAXIT ) );
            }
            S1 = ONE;
            WR = ESHIFT;

         } else {

            // Shifts based on the generalized eigenvalues of the
            // bottom-right 2x2 block of A and B. The first eigenvalue
            // returned by SLAG2 is the Wilkinson shift (AEP p.512),

            slag2(H( ILAST-1, ILAST-1 ), LDH, T( ILAST-1, ILAST-1 ), LDT, SAFMIN*SAFETY, S1, S2, WR, WR2, WI );

            if ( ABS( (WR/S1)*T( ILAST, ILAST ) - H( ILAST, ILAST ) ) > ABS( (WR2/S2)*T( ILAST, ILAST ) - H( ILAST, ILAST ) ) ) {
               TEMP = WR;
               WR = WR2;
               WR2 = TEMP;
               TEMP = S1;
               S1 = S2;
               S2 = TEMP;
            }
            TEMP = max( S1, SAFMIN*max( ONE, ( WR ).abs(), ( WI ).abs() ) );
            if (WI != ZERO) GO TO 200;
         }

         // Fiddle with shift to avoid overflow

         TEMP = min( ASCALE, ONE )*( HALF*SAFMAX );
         if ( S1 > TEMP ) {
            SCALE = TEMP / S1;
         } else {
            SCALE = ONE;
         }

         TEMP = min( BSCALE, ONE )*( HALF*SAFMAX );
         if( ( WR ).abs() > TEMP ) SCALE = min( SCALE, TEMP / ( WR ).abs() );
         S1 = SCALE*S1;
         WR = SCALE*WR;

         // Now check for two consecutive small subdiagonals.

         for (J = ILAST - 1; J >= IFIRST + 1; J--) { // 120
            ISTART = J;
            TEMP = ( S1*H( J, J-1 ) ).abs();
            TEMP2 = ABS( S1*H( J, J )-WR*T( J, J ) );
            TEMPR = max( TEMP, TEMP2 );
            if ( TEMPR < ONE && TEMPR != ZERO ) {
               TEMP = TEMP / TEMPR;
               TEMP2 = TEMP2 / TEMPR;
            }
            if( ABS( ( ASCALE*H( J+1, J ) )*TEMP ) <= ( ASCALE*ATOL )* TEMP2 )GO TO 130;
         } // 120

         ISTART = IFIRST;
         } // 130

         // Do an implicit single-shift QZ sweep.

         // Initial Q

         TEMP = S1*H( ISTART, ISTART ) - WR*T( ISTART, ISTART );
         TEMP2 = S1*H( ISTART+1, ISTART );
         slartg(TEMP, TEMP2, C, S, TEMPR );

         // Sweep

         for (J = ISTART; J <= ILAST - 1; J++) { // 190
            if ( J > ISTART ) {
               TEMP = H( J, J-1 );
               slartg(TEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) );
               H[J+1][J-1] = ZERO;
            }

            for (JC = J; JC <= ILASTM; JC++) { // 140
               TEMP = C*H( J, JC ) + S*H( J+1, JC );
               H[J+1][JC] = -S*H( J, JC ) + C*H( J+1, JC );
               H[J][JC] = TEMP;
               TEMP2 = C*T( J, JC ) + S*T( J+1, JC );
               T[J+1][JC] = -S*T( J, JC ) + C*T( J+1, JC );
               T[J][JC] = TEMP2;
            } // 140
            if ( ILQ ) {
               for (JR = 1; JR <= N; JR++) { // 150
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 );
                  Q[JR][J+1] = -S*Q( JR, J ) + C*Q( JR, J+1 );
                  Q[JR][J] = TEMP;
               } // 150
            }

            TEMP = T( J+1, J+1 );
            slartg(TEMP, T( J+1, J ), C, S, T( J+1, J+1 ) );
            T[J+1][J] = ZERO;

            for (JR = IFRSTM; JR <= min( J+2, ILAST ); JR++) { // 160
               TEMP = C*H( JR, J+1 ) + S*H( JR, J );
               H[JR][J] = -S*H( JR, J+1 ) + C*H( JR, J );
               H[JR][J+1] = TEMP;
            } // 160
            for (JR = IFRSTM; JR <= J; JR++) { // 170
               TEMP = C*T( JR, J+1 ) + S*T( JR, J );
               T[JR][J] = -S*T( JR, J+1 ) + C*T( JR, J );
               T[JR][J+1] = TEMP;
            } // 170
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 180
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J );
                  Z[JR][J] = -S*Z( JR, J+1 ) + C*Z( JR, J );
                  Z[JR][J+1] = TEMP;
               } // 180
            }
         } // 190

         GO TO 350;

         // Use Francis double-shift

         // Note: the Francis double-shift should work with real shifts,
         //       but only if the block is at least 3x3.
         //       This code may break if this point is reached with
         //       a 2x2 block with real eigenvalues.

         } // 200
         if ( IFIRST+1 == ILAST ) {

            // Special case -- 2x2 block with complex eigenvectors

            // Step 1: Standardize, that is, rotate so that

                        // ( B11  0  )
                    // B = (         )  with B11 non-negative.
                    //     (  0  B22 )

            slasv2(T( ILAST-1, ILAST-1 ), T( ILAST-1, ILAST ), T( ILAST, ILAST ), B22, B11, SR, CR, SL, CL );

            if ( B11 < ZERO ) {
               CR = -CR;
               SR = -SR;
               B11 = -B11;
               B22 = -B22;
            }

            srot(ILASTM+1-IFIRST, H( ILAST-1, ILAST-1 ), LDH, H( ILAST, ILAST-1 ), LDH, CL, SL );
            srot(ILAST+1-IFRSTM, H( IFRSTM, ILAST-1 ), 1, H( IFRSTM, ILAST ), 1, CR, SR );

            if (ILAST < ILASTM) srot( ILASTM-ILAST, T( ILAST-1, ILAST+1 ), LDT, T( ILAST, ILAST+1 ), LDT, CL, SL );
            IF( IFRSTM < ILAST-1 ) srot( IFIRST-IFRSTM, T( IFRSTM, ILAST-1 ), 1, T( IFRSTM, ILAST ), 1, CR, SR );

            if (ILQ) srot( N, Q( 1, ILAST-1 ), 1, Q( 1, ILAST ), 1, CL, SL );
            IF( ILZ ) srot( N, Z( 1, ILAST-1 ), 1, Z( 1, ILAST ), 1, CR, SR );

            T[ILAST-1][ILAST-1] = B11;
            T[ILAST-1][ILAST] = ZERO;
            T[ILAST][ILAST-1] = ZERO;
            T[ILAST][ILAST] = B22;

            // If B22 is negative, negate column ILAST

            if ( B22 < ZERO ) {
               for (J = IFRSTM; J <= ILAST; J++) { // 210
                  H[J][ILAST] = -H( J, ILAST );
                  T[J][ILAST] = -T( J, ILAST );
               } // 210

               if ( ILZ ) {
                  for (J = 1; J <= N; J++) { // 220
                     Z[J][ILAST] = -Z( J, ILAST );
                  } // 220
               }
               B22 = -B22;
            }

            // Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)

            // Recompute shift

            slag2(H( ILAST-1, ILAST-1 ), LDH, T( ILAST-1, ILAST-1 ), LDT, SAFMIN*SAFETY, S1, TEMP, WR, TEMP2, WI );

            // If standardization has perturbed the shift onto real line,
            // do another (real single-shift) QR step.

            if (WI == ZERO) GO TO 350;
            S1INV = ONE / S1;

            // Do EISPACK (QZVAL) computation of alpha and beta

            A11 = H( ILAST-1, ILAST-1 );
            A21 = H( ILAST, ILAST-1 );
            A12 = H( ILAST-1, ILAST );
            A22 = H( ILAST, ILAST );

            // Compute complex Givens rotation on right
            // (Assume some element of C = (sA - wB) > unfl )
            //                  __
            // (sA - wB) ( CZ   -SZ )
            //           ( SZ    CZ )

            C11R = S1*A11 - WR*B11;
            C11I = -WI*B11;
            C12 = S1*A12;
            C21 = S1*A21;
            C22R = S1*A22 - WR*B22;
            C22I = -WI*B22;

            if ( ( C11R ).abs()+( C11I ).abs()+( C12 ).abs() > ( C21 ).abs()+ ( C22R ).abs()+( C22I ).abs() ) {
               T1 = SLAPY3( C12, C11R, C11I );
               CZ = C12 / T1;
               SZR = -C11R / T1;
               SZI = -C11I / T1;
            } else {
               CZ = SLAPY2( C22R, C22I );
               if ( CZ <= SAFMIN ) {
                  CZ = ZERO;
                  SZR = ONE;
                  SZI = ZERO;
               } else {
                  TEMPR = C22R / CZ;
                  TEMPI = C22I / CZ;
                  T1 = SLAPY2( CZ, C21 );
                  CZ = CZ / T1;
                  SZR = -C21*TEMPR / T1;
                  SZI = C21*TEMPI / T1;
               }
            }

            // Compute Givens rotation on left

            // (  CQ   SQ )
            // (  __      )  A or B
            // ( -SQ   CQ )

            AN = ( A11 ).abs() + ( A12 ).abs() + ( A21 ).abs() + ( A22 ).abs();
            BN = ( B11 ).abs() + ( B22 ).abs();
            WABS = ( WR ).abs() + ( WI ).abs();
            if ( S1*AN > WABS*BN ) {
               CQ = CZ*B11;
               SQR = SZR*B22;
               SQI = -SZI*B22;
            } else {
               A1R = CZ*A11 + SZR*A12;
               A1I = SZI*A12;
               A2R = CZ*A21 + SZR*A22;
               A2I = SZI*A22;
               CQ = SLAPY2( A1R, A1I );
               if ( CQ <= SAFMIN ) {
                  CQ = ZERO;
                  SQR = ONE;
                  SQI = ZERO;
               } else {
                  TEMPR = A1R / CQ;
                  TEMPI = A1I / CQ;
                  SQR = TEMPR*A2R + TEMPI*A2I;
                  SQI = TEMPI*A2R - TEMPR*A2I;
               }
            }
            T1 = SLAPY3( CQ, SQR, SQI );
            CQ = CQ / T1;
            SQR = SQR / T1;
            SQI = SQI / T1;

            // Compute diagonal elements of QBZ

            TEMPR = SQR*SZR - SQI*SZI;
            TEMPI = SQR*SZI + SQI*SZR;
            B1R = CQ*CZ*B11 + TEMPR*B22;
            B1I = TEMPI*B22;
            B1A = SLAPY2( B1R, B1I );
            B2R = CQ*CZ*B22 + TEMPR*B11;
            B2I = -TEMPI*B11;
            B2A = SLAPY2( B2R, B2I );

            // Normalize so beta > 0, and Im( alpha1 ) > 0

            BETA[ILAST-1] = B1A;
            BETA[ILAST] = B2A;
            ALPHAR[ILAST-1] = ( WR*B1A )*S1INV;
            ALPHAI[ILAST-1] = ( WI*B1A )*S1INV;
            ALPHAR[ILAST] = ( WR*B2A )*S1INV;
            ALPHAI[ILAST] = -( WI*B2A )*S1INV;

            // Step 3: Go to next block -- exit if finished.

            ILAST = IFIRST - 1;
            if (ILAST < ILO) GO TO 380;

            // Reset counters

            IITER = 0;
            ESHIFT = ZERO;
            if ( !ILSCHR ) {
               ILASTM = ILAST;
               if (IFRSTM > ILAST) IFRSTM = ILO;
            }
            GO TO 350;
         } else {

            // Usual case: 3x3 or larger block, using Francis implicit
            //             double-shift

                                     // 2
            // Eigenvalue equation is  w  - c w + d = 0,

                                          // -1 2        -1
            // so compute 1st column of  (A B  )  - c A B   + d
            // using the formula in QZIT (from EISPACK)

            // We assume that the block is at least 3x3

            AD11 = ( ASCALE*H( ILAST-1, ILAST-1 ) ) / ( BSCALE*T( ILAST-1, ILAST-1 ) )             AD21 = ( ASCALE*H( ILAST, ILAST-1 ) ) / ( BSCALE*T( ILAST-1, ILAST-1 ) )             AD12 = ( ASCALE*H( ILAST-1, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) )             AD22 = ( ASCALE*H( ILAST, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) );
            U12 = T( ILAST-1, ILAST ) / T( ILAST, ILAST );
            AD11L = ( ASCALE*H( IFIRST, IFIRST ) ) / ( BSCALE*T( IFIRST, IFIRST ) )             AD21L = ( ASCALE*H( IFIRST+1, IFIRST ) ) / ( BSCALE*T( IFIRST, IFIRST ) )             AD12L = ( ASCALE*H( IFIRST, IFIRST+1 ) ) / ( BSCALE*T( IFIRST+1, IFIRST+1 ) )             AD22L = ( ASCALE*H( IFIRST+1, IFIRST+1 ) ) / ( BSCALE*T( IFIRST+1, IFIRST+1 ) )             AD32L = ( ASCALE*H( IFIRST+2, IFIRST+1 ) ) / ( BSCALE*T( IFIRST+1, IFIRST+1 ) );
            U12L = T( IFIRST, IFIRST+1 ) / T( IFIRST+1, IFIRST+1 );

            V[1] = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 + AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L             V( 2 ) = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )- ( AD22-AD11L )+AD21*U12 )*AD21L;
            V[3] = AD32L*AD21L;

            ISTART = IFIRST;

            slarfg(3, V( 1 ), V( 2 ), 1, TAU );
            V[1] = ONE;

            // Sweep

            for (J = ISTART; J <= ILAST - 2; J++) { // 290

               // All but last elements: use 3x3 Householder transforms.

               // Zero (j-1)st column of A

               if ( J > ISTART ) {
                  V[1] = H( J, J-1 );
                  V[2] = H( J+1, J-1 );
                  V[3] = H( J+2, J-1 );

                  slarfg(3, H( J, J-1 ), V( 2 ), 1, TAU );
                  V[1] = ONE;
                  H[J+1][J-1] = ZERO;
                  H[J+2][J-1] = ZERO;
               }

               T2 = TAU * V( 2 );
               T3 = TAU * V( 3 );
               for (JC = J; JC <= ILASTM; JC++) { // 230
                  TEMP = H( J, JC )+V( 2 )*H( J+1, JC )+V( 3 )* H( J+2, JC );
                  H[J][JC] = H( J, JC ) - TEMP*TAU;
                  H[J+1][JC] = H( J+1, JC ) - TEMP*T2;
                  H[J+2][JC] = H( J+2, JC ) - TEMP*T3;
                  TEMP2 = T( J, JC )+V( 2 )*T( J+1, JC )+V( 3 )* T( J+2, JC );
                  T[J][JC] = T( J, JC ) - TEMP2*TAU;
                  T[J+1][JC] = T( J+1, JC ) - TEMP2*T2;
                  T[J+2][JC] = T( J+2, JC ) - TEMP2*T3;
               } // 230
               if ( ILQ ) {
                  for (JR = 1; JR <= N; JR++) { // 240
                     TEMP = Q( JR, J )+V( 2 )*Q( JR, J+1 )+V( 3 )* Q( JR, J+2 );
                     Q[JR][J] = Q( JR, J ) - TEMP*TAU;
                     Q[JR][J+1] = Q( JR, J+1 ) - TEMP*T2;
                     Q[JR][J+2] = Q( JR, J+2 ) - TEMP*T3;
                  } // 240
               }

               // Zero j-th column of B (see SLAGBC for details)

               // Swap rows to pivot

               ILPIVT = false;
               TEMP = max( ( T( J+1, J+1 ) ).abs(), ( T( J+1, J+2 ) ).abs() );
               TEMP2 = max( ( T( J+2, J+1 ) ).abs(), ( T( J+2, J+2 ) ).abs() );
               if ( max( TEMP, TEMP2 ) < SAFMIN ) {
                  SCALE = ZERO;
                  U1 = ONE;
                  U2 = ZERO;
                  GO TO 250;
               } else if ( TEMP >= TEMP2 ) {
                  W11 = T( J+1, J+1 );
                  W21 = T( J+2, J+1 );
                  W12 = T( J+1, J+2 );
                  W22 = T( J+2, J+2 );
                  U1 = T( J+1, J );
                  U2 = T( J+2, J );
               } else {
                  W21 = T( J+1, J+1 );
                  W11 = T( J+2, J+1 );
                  W22 = T( J+1, J+2 );
                  W12 = T( J+2, J+2 );
                  U2 = T( J+1, J );
                  U1 = T( J+2, J );
               }

               // Swap columns if nec.

               if ( ( W12 ).abs() > ( W11 ).abs() ) {
                  ILPIVT = true;
                  TEMP = W12;
                  TEMP2 = W22;
                  W12 = W11;
                  W22 = W21;
                  W11 = TEMP;
                  W21 = TEMP2;
               }

               // LU-factor

               TEMP = W21 / W11;
               U2 = U2 - TEMP*U1;
               W22 = W22 - TEMP*W12;
               W21 = ZERO;

               // Compute SCALE

               SCALE = ONE;
               if ( ( W22 ).abs() < SAFMIN ) {
                  SCALE = ZERO;
                  U2 = ONE;
                  U1 = -W12 / W11;
                  GO TO 250;
               }
               if( ( W22 ).abs() < ( U2 ).abs() ) SCALE = ( W22 / U2 ).abs();
               IF( ( W11 ).abs() < ( U1 ).abs() ) SCALE = min( SCALE, ( W11 / U1 ).abs() );

               // Solve

               U2 = ( SCALE*U2 ) / W22;
               U1 = ( SCALE*U1-W12*U2 ) / W11;

               } // 250
               if ( ILPIVT ) {
                  TEMP = U2;
                  U2 = U1;
                  U1 = TEMP;
               }

               // Compute Householder Vector

               T1 = sqrt( SCALE**2+U1**2+U2**2 );
               TAU = ONE + SCALE / T1;
               VS = -ONE / ( SCALE+T1 );
               V[1] = ONE;
               V[2] = VS*U1;
               V[3] = VS*U2;

               // Apply transformations from the right.

               T2 = TAU*V( 2 );
               T3 = TAU*V( 3 );
               for (JR = IFRSTM; JR <= min( J+3, ILAST ); JR++) { // 260
                  TEMP = H( JR, J )+V( 2 )*H( JR, J+1 )+V( 3 )* H( JR, J+2 );
                  H[JR][J] = H( JR, J ) - TEMP*TAU;
                  H[JR][J+1] = H( JR, J+1 ) - TEMP*T2;
                  H[JR][J+2] = H( JR, J+2 ) - TEMP*T3;
               } // 260
               for (JR = IFRSTM; JR <= J + 2; JR++) { // 270
                  TEMP = T( JR, J )+V( 2 )*T( JR, J+1 )+V( 3 )* T( JR, J+2 );
                  T[JR][J] = T( JR, J ) - TEMP*TAU;
                  T[JR][J+1] = T( JR, J+1 ) - TEMP*T2;
                  T[JR][J+2] = T( JR, J+2 ) - TEMP*T3;
               } // 270
               if ( ILZ ) {
                  for (JR = 1; JR <= N; JR++) { // 280
                     TEMP = Z( JR, J )+V( 2 )*Z( JR, J+1 )+V( 3 )* Z( JR, J+2 );
                     Z[JR][J] = Z( JR, J ) - TEMP*TAU;
                     Z[JR][J+1] = Z( JR, J+1 ) - TEMP*T2;
                     Z[JR][J+2] = Z( JR, J+2 ) - TEMP*T3;
                  } // 280
               }
               T[J+1][J] = ZERO;
               T[J+2][J] = ZERO;
            } // 290

            // Last elements: Use Givens rotations

            // Rotations from the left

            J = ILAST - 1;
            TEMP = H( J, J-1 );
            slartg(TEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) );
            H[J+1][J-1] = ZERO;

            for (JC = J; JC <= ILASTM; JC++) { // 300
               TEMP = C*H( J, JC ) + S*H( J+1, JC );
               H[J+1][JC] = -S*H( J, JC ) + C*H( J+1, JC );
               H[J][JC] = TEMP;
               TEMP2 = C*T( J, JC ) + S*T( J+1, JC );
               T[J+1][JC] = -S*T( J, JC ) + C*T( J+1, JC );
               T[J][JC] = TEMP2;
            } // 300
            if ( ILQ ) {
               for (JR = 1; JR <= N; JR++) { // 310
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 );
                  Q[JR][J+1] = -S*Q( JR, J ) + C*Q( JR, J+1 );
                  Q[JR][J] = TEMP;
               } // 310
            }

            // Rotations from the right.

            TEMP = T( J+1, J+1 );
            slartg(TEMP, T( J+1, J ), C, S, T( J+1, J+1 ) );
            T[J+1][J] = ZERO;

            for (JR = IFRSTM; JR <= ILAST; JR++) { // 320
               TEMP = C*H( JR, J+1 ) + S*H( JR, J );
               H[JR][J] = -S*H( JR, J+1 ) + C*H( JR, J );
               H[JR][J+1] = TEMP;
            } // 320
            for (JR = IFRSTM; JR <= ILAST - 1; JR++) { // 330
               TEMP = C*T( JR, J+1 ) + S*T( JR, J );
               T[JR][J] = -S*T( JR, J+1 ) + C*T( JR, J );
               T[JR][J+1] = TEMP;
            } // 330
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 340
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J );
                  Z[JR][J] = -S*Z( JR, J+1 ) + C*Z( JR, J );
                  Z[JR][J+1] = TEMP;
               } // 340
            }

            // End of Double-Shift code

         }

         GO TO 350;

         // End of iteration loop

         } // 350
      } // 360

      // Drop-through = non-convergence

      INFO = ILAST;
      GO TO 420;

      // Successful completion of all QZ steps

      } // 380

      // Set Eigenvalues 1:ILO-1

      for (J = 1; J <= ILO - 1; J++) { // 410
         if ( T( J, J ) < ZERO ) {
            if ( ILSCHR ) {
               for (JR = 1; JR <= J; JR++) { // 390
                  H[JR][J] = -H( JR, J );
                  T[JR][J] = -T( JR, J );
               } // 390
            } else {
               H[J][J] = -H( J, J );
               T[J][J] = -T( J, J );
            }
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 400
                  Z[JR][J] = -Z( JR, J );
               } // 400
            }
         }
         ALPHAR[J] = H( J, J );
         ALPHAI[J] = ZERO;
         BETA[J] = T( J, J );
      } // 410

      // Normal Termination

      INFO = 0;

      // Exit (other than argument error) -- return optimal workspace size

      } // 420
      WORK[1] = SROUNDUP_LWORK( N );
      }

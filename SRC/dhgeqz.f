      SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, COMPZ, JOB;
      int                IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             ALPHAI( * ), ALPHAR( * ), BETA( * ), H( LDH, * ), Q( LDQ, * ), T( LDT, * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
*    $                     SAFETY = 1.0E+0 )
      double             HALF, ZERO, ONE, SAFETY;
      const              HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0, SAFETY = 1.0D+2 ;
      // ..
      // .. Local Scalars ..
      bool               ILAZR2, ILAZRO, ILPIVT, ILQ, ILSCHR, ILZ, LQUERY;
      int                ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST, ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER, JR, MAXIT;
      double             A11, A12, A1I, A1R, A21, A22, A2I, A2R, AD11, AD11L, AD12, AD12L, AD21, AD21L, AD22, AD22L, AD32L, AN, ANORM, ASCALE, ATOL, B11, B1A, B1I, B1R, B22, B2A, B2I, B2R, BN, BNORM, BSCALE, BTOL, C, C11I, C11R, C12, C21, C22I, C22R, CL, CQ, CR, CZ, ESHIFT, S, S1, S1INV, S2, SAFMAX, SAFMIN, SCALE, SL, SQI, SQR, SR, SZI, SZR, T1, T2, T3, TAU, TEMP, TEMP2, TEMPI, TEMPR, U1, U12, U12L, U2, ULP, VS, W11, W12, W21, W22, WABS, WI, WR, WR2;
      // ..
      // .. Local Arrays ..
      double             V( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANHS, DLAPY2, DLAPY3;
      // EXTERNAL LSAME, DLAMCH, DLANHS, DLAPY2, DLAPY3
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAG2, DLARFG, DLARTG, DLASET, DLASV2, DROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Decode JOB, COMPQ, COMPZ

      if ( LSAME( JOB, 'E' ) ) {
         ILSCHR = .FALSE.
         ISCHUR = 1
      } else if ( LSAME( JOB, 'S' ) ) {
         ILSCHR = .TRUE.
         ISCHUR = 2
      } else {
         ISCHUR = 0
      }

      if ( LSAME( COMPQ, 'N' ) ) {
         ILQ = .FALSE.
         ICOMPQ = 1
      } else if ( LSAME( COMPQ, 'V' ) ) {
         ILQ = .TRUE.
         ICOMPQ = 2
      } else if ( LSAME( COMPQ, 'I' ) ) {
         ILQ = .TRUE.
         ICOMPQ = 3
      } else {
         ICOMPQ = 0
      }

      if ( LSAME( COMPZ, 'N' ) ) {
         ILZ = .FALSE.
         ICOMPZ = 1
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ILZ = .TRUE.
         ICOMPZ = 2
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ILZ = .TRUE.
         ICOMPZ = 3
      } else {
         ICOMPZ = 0
      }

      // Check Argument Values

      INFO = 0
      WORK( 1 ) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 )
      if ( ISCHUR.EQ.0 ) {
         INFO = -1
      } else if ( ICOMPQ.EQ.0 ) {
         INFO = -2
      } else if ( ICOMPZ.EQ.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( ILO.LT.1 ) {
         INFO = -5
      } else if ( IHI.GT.N .OR. IHI.LT.ILO-1 ) {
         INFO = -6
      } else if ( LDH.LT.N ) {
         INFO = -8
      } else if ( LDT.LT.N ) {
         INFO = -10
      } else if ( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) {
         INFO = -15
      } else if ( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) {
         INFO = -17
      } else if ( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
         INFO = -19
      }
      if ( INFO.NE.0 ) {
         xerbla('DHGEQZ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.LE.0 ) {
         WORK( 1 ) = DBLE( 1 )
         RETURN
      }

      // Initialize Q and Z

      IF( ICOMPQ.EQ.3 ) CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )       IF( ICOMPZ.EQ.3 ) CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )

      // Machine Constants

      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = DLANHS( 'F', IN, H( ILO, ILO ), LDH, WORK )
      BNORM = DLANHS( 'F', IN, T( ILO, ILO ), LDT, WORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )

      // Set Eigenvalues IHI+1:N

      for (J = IHI + 1; J <= N; J++) { // 30
         if ( T( J, J ).LT.ZERO ) {
            if ( ILSCHR ) {
               for (JR = 1; JR <= J; JR++) { // 10
                  H( JR, J ) = -H( JR, J )
                  T( JR, J ) = -T( JR, J )
               } // 10
            } else {
               H( J, J ) = -H( J, J )
               T( J, J ) = -T( J, J )
            }
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 20
                  Z( JR, J ) = -Z( JR, J )
               } // 20
            }
         }
         ALPHAR( J ) = H( J, J )
         ALPHAI( J ) = ZERO
         BETA( J ) = T( J, J )
      } // 30

      // If IHI < ILO, skip QZ steps

      IF( IHI.LT.ILO ) GO TO 380

      // MAIN QZ ITERATION LOOP

      // Initialize dynamic indices

      // Eigenvalues ILAST+1:N have been found.
         // Column operations modify rows IFRSTM:whatever.
         // Row operations modify columns whatever:ILASTM.

      // If only eigenvalues are being computed, then
         // IFRSTM is the row of the last splitting row above row ILAST;
         // this is always at least ILO.
      // IITER counts iterations since the last eigenvalue was found,
         // to tell when to use an extraordinary shift.
      // MAXIT is the maximum number of QZ sweeps allowed.

      ILAST = IHI
      if ( ILSCHR ) {
         IFRSTM = 1
         ILASTM = N
      } else {
         IFRSTM = ILO
         ILASTM = IHI
      }
      IITER = 0
      ESHIFT = ZERO
      MAXIT = 30*( IHI-ILO+1 )

      for (JITER = 1; JITER <= MAXIT; JITER++) { // 360

         // Split the matrix if possible.

         // Two tests:
            // 1: H(j,j-1)=0  or  j=ILO
            // 2: T(j,j)=0

         if ( ILAST.EQ.ILO ) {

            // Special case: j=ILAST

            GO TO 80
         } else {
            if ( ABS( H( ILAST, ILAST-1 ) ).LE.MAX( SAFMIN, ULP*(  ABS( H( ILAST, ILAST ) ) + ABS( H( ILAST-1, ILAST-1 ) ) ) ) ) {
               H( ILAST, ILAST-1 ) = ZERO
               GO TO 80
            }
         }

         if ( ABS( T( ILAST, ILAST ) ).LE.BTOL ) {
            T( ILAST, ILAST ) = ZERO
            GO TO 70
         }

         // General case: j<ILAST

         DO 60 J = ILAST - 1, ILO, -1

            // Test 1: for H(j,j-1)=0 or j=ILO

            if ( J.EQ.ILO ) {
               ILAZRO = .TRUE.
            } else {
               if ( ABS( H( J, J-1 ) ).LE.MAX( SAFMIN, ULP*(  ABS( H( J, J ) ) + ABS( H( J-1, J-1 ) ) ) ) ) {
                  H( J, J-1 ) = ZERO
                  ILAZRO = .TRUE.
               } else {
                  ILAZRO = .FALSE.
               }
            }

            // Test 2: for T(j,j)=0

            if ( ABS( T( J, J ) ).LT.BTOL ) {
               T( J, J ) = ZERO

               // Test 1a: Check for 2 consecutive small subdiagonals in A

               ILAZR2 = .FALSE.
               if ( .NOT.ILAZRO ) {
                  TEMP = ABS( H( J, J-1 ) )
                  TEMP2 = ABS( H( J, J ) )
                  TEMPR = MAX( TEMP, TEMP2 )
                  if ( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) {
                     TEMP = TEMP / TEMPR
                     TEMP2 = TEMP2 / TEMPR
                  }
                  IF( TEMP*( ASCALE*ABS( H( J+1, J ) ) ).LE.TEMP2* ( ASCALE*ATOL ) )ILAZR2 = .TRUE.
               }

               // If both tests pass (1 & 2), i.e., the leading diagonal
               // element of B in the block is zero, split a 1x1 block off
               // at the top. (I.e., at the J-th row/column) The leading
               // diagonal element of the remainder can also be zero, so
               // this may have to be done repeatedly.

               if ( ILAZRO .OR. ILAZR2 ) {
                  for (JCH = J; JCH <= ILAST - 1; JCH++) { // 40
                     TEMP = H( JCH, JCH )
                     dlartg(TEMP, H( JCH+1, JCH ), C, S, H( JCH, JCH ) );
                     H( JCH+1, JCH ) = ZERO
                     drot(ILASTM-JCH, H( JCH, JCH+1 ), LDH, H( JCH+1, JCH+1 ), LDH, C, S )                      CALL DROT( ILASTM-JCH, T( JCH, JCH+1 ), LDT, T( JCH+1, JCH+1 ), LDT, C, S )                      IF( ILQ ) CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1, C, S );
                     IF( ILAZR2 ) H( JCH, JCH-1 ) = H( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     if ( ABS( T( JCH+1, JCH+1 ) ).GE.BTOL ) {
                        if ( JCH+1.GE.ILAST ) {
                           GO TO 80
                        } else {
                           IFIRST = JCH + 1
                           GO TO 110
                        }
                     }
                     T( JCH+1, JCH+1 ) = ZERO
                  } // 40
                  GO TO 70
               } else {

                  // Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                  // Then process as in the case T(ILAST,ILAST)=0

                  for (JCH = J; JCH <= ILAST - 1; JCH++) { // 50
                     TEMP = T( JCH, JCH+1 )
                     dlartg(TEMP, T( JCH+1, JCH+1 ), C, S, T( JCH, JCH+1 ) );
                     T( JCH+1, JCH+1 ) = ZERO
                     IF( JCH.LT.ILASTM-1 ) CALL DROT( ILASTM-JCH-1, T( JCH, JCH+2 ), LDT, T( JCH+1, JCH+2 ), LDT, C, S )                      CALL DROT( ILASTM-JCH+2, H( JCH, JCH-1 ), LDH, H( JCH+1, JCH-1 ), LDH, C, S )                      IF( ILQ ) CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1, C, S )
                     TEMP = H( JCH+1, JCH )
                     dlartg(TEMP, H( JCH+1, JCH-1 ), C, S, H( JCH+1, JCH ) );
                     H( JCH+1, JCH-1 ) = ZERO
                     drot(JCH+1-IFRSTM, H( IFRSTM, JCH ), 1, H( IFRSTM, JCH-1 ), 1, C, S )                      CALL DROT( JCH-IFRSTM, T( IFRSTM, JCH ), 1, T( IFRSTM, JCH-1 ), 1, C, S )                      IF( ILZ ) CALL DROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1, C, S );
                  } // 50
                  GO TO 70
               }
            } else if ( ILAZRO ) {

               // Only test 1 passed -- work on J:ILAST

               IFIRST = J
               GO TO 110
            }

            // Neither test passed -- try next J

         } // 60

         // (Drop-through is "impossible")

         INFO = N + 1
         GO TO 420

         // T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
         // 1x1 block.

         } // 70
         TEMP = H( ILAST, ILAST )
         dlartg(TEMP, H( ILAST, ILAST-1 ), C, S, H( ILAST, ILAST ) );
         H( ILAST, ILAST-1 ) = ZERO
         drot(ILAST-IFRSTM, H( IFRSTM, ILAST ), 1, H( IFRSTM, ILAST-1 ), 1, C, S )          CALL DROT( ILAST-IFRSTM, T( IFRSTM, ILAST ), 1, T( IFRSTM, ILAST-1 ), 1, C, S )          IF( ILZ ) CALL DROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S );

         // H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
                               // and BETA

         } // 80
         if ( T( ILAST, ILAST ).LT.ZERO ) {
            if ( ILSCHR ) {
               for (J = IFRSTM; J <= ILAST; J++) { // 90
                  H( J, ILAST ) = -H( J, ILAST )
                  T( J, ILAST ) = -T( J, ILAST )
               } // 90
            } else {
               H( ILAST, ILAST ) = -H( ILAST, ILAST )
               T( ILAST, ILAST ) = -T( ILAST, ILAST )
            }
            if ( ILZ ) {
               for (J = 1; J <= N; J++) { // 100
                  Z( J, ILAST ) = -Z( J, ILAST )
               } // 100
            }
         }
         ALPHAR( ILAST ) = H( ILAST, ILAST )
         ALPHAI( ILAST ) = ZERO
         BETA( ILAST ) = T( ILAST, ILAST )

         // Go to next block -- exit if finished.

         ILAST = ILAST - 1
         IF( ILAST.LT.ILO ) GO TO 380

         // Reset counters

         IITER = 0
         ESHIFT = ZERO
         if ( .NOT.ILSCHR ) {
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST ) IFRSTM = ILO
         }
         GO TO 350

         // QZ step

         // This iteration only involves rows/columns IFIRST:ILAST. We
         // assume IFIRST < ILAST, and that the diagonal of B is non-zero.

         } // 110
         IITER = IITER + 1
         if ( .NOT.ILSCHR ) {
            IFRSTM = IFIRST
         }

         // Compute single shifts.

         // At this point, IFIRST < ILAST, and the diagonal elements of
         // T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
         // magnitude)

         if ( ( IITER / 10 )*10.EQ.IITER ) {

            // Exceptional shift.  Chosen for no particularly good reason.
            // (Single shift only.)

            IF( ( DBLE( MAXIT )*SAFMIN )*ABS( H( ILAST, ILAST-1 ) ).LT. ABS( T( ILAST-1, ILAST-1 ) ) ) THEN                ESHIFT = H( ILAST, ILAST-1 ) / T( ILAST-1, ILAST-1 )
            } else {
               ESHIFT = ESHIFT + ONE / ( SAFMIN*DBLE( MAXIT ) )
            }
            S1 = ONE
            WR = ESHIFT

         } else {

            // Shifts based on the generalized eigenvalues of the
            // bottom-right 2x2 block of A and B. The first eigenvalue
            // returned by DLAG2 is the Wilkinson shift (AEP p.512),

            dlag2(H( ILAST-1, ILAST-1 ), LDH, T( ILAST-1, ILAST-1 ), LDT, SAFMIN*SAFETY, S1, S2, WR, WR2, WI );

            if ( ABS( (WR/S1)*T( ILAST, ILAST ) - H( ILAST, ILAST ) ) .GT. ABS( (WR2/S2)*T( ILAST, ILAST ) - H( ILAST, ILAST ) ) ) {
               TEMP = WR
               WR = WR2
               WR2 = TEMP
               TEMP = S1
               S1 = S2
               S2 = TEMP
            }
            TEMP = MAX( S1, SAFMIN*MAX( ONE, ABS( WR ), ABS( WI ) ) )
            IF( WI.NE.ZERO ) GO TO 200
         }

         // Fiddle with shift to avoid overflow

         TEMP = MIN( ASCALE, ONE )*( HALF*SAFMAX )
         if ( S1.GT.TEMP ) {
            SCALE = TEMP / S1
         } else {
            SCALE = ONE
         }

         TEMP = MIN( BSCALE, ONE )*( HALF*SAFMAX )
         IF( ABS( WR ).GT.TEMP ) SCALE = MIN( SCALE, TEMP / ABS( WR ) )
         S1 = SCALE*S1
         WR = SCALE*WR

         // Now check for two consecutive small subdiagonals.

         DO 120 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            TEMP = ABS( S1*H( J, J-1 ) )
            TEMP2 = ABS( S1*H( J, J )-WR*T( J, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            if ( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) {
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            }
            IF( ABS( ( ASCALE*H( J+1, J ) )*TEMP ).LE.( ASCALE*ATOL )* TEMP2 )GO TO 130
         } // 120

         ISTART = IFIRST
         } // 130

         // Do an implicit single-shift QZ sweep.

         // Initial Q

         TEMP = S1*H( ISTART, ISTART ) - WR*T( ISTART, ISTART )
         TEMP2 = S1*H( ISTART+1, ISTART )
         dlartg(TEMP, TEMP2, C, S, TEMPR );

         // Sweep

         for (J = ISTART; J <= ILAST - 1; J++) { // 190
            if ( J.GT.ISTART ) {
               TEMP = H( J, J-1 )
               dlartg(TEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) );
               H( J+1, J-1 ) = ZERO
            }

            for (JC = J; JC <= ILASTM; JC++) { // 140
               TEMP = C*H( J, JC ) + S*H( J+1, JC )
               H( J+1, JC ) = -S*H( J, JC ) + C*H( J+1, JC )
               H( J, JC ) = TEMP
               TEMP2 = C*T( J, JC ) + S*T( J+1, JC )
               T( J+1, JC ) = -S*T( J, JC ) + C*T( J+1, JC )
               T( J, JC ) = TEMP2
            } // 140
            if ( ILQ ) {
               for (JR = 1; JR <= N; JR++) { // 150
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = TEMP
               } // 150
            }

            TEMP = T( J+1, J+1 )
            dlartg(TEMP, T( J+1, J ), C, S, T( J+1, J+1 ) );
            T( J+1, J ) = ZERO

            DO 160 JR = IFRSTM, MIN( J+2, ILAST )
               TEMP = C*H( JR, J+1 ) + S*H( JR, J )
               H( JR, J ) = -S*H( JR, J+1 ) + C*H( JR, J )
               H( JR, J+1 ) = TEMP
            } // 160
            for (JR = IFRSTM; JR <= J; JR++) { // 170
               TEMP = C*T( JR, J+1 ) + S*T( JR, J )
               T( JR, J ) = -S*T( JR, J+1 ) + C*T( JR, J )
               T( JR, J+1 ) = TEMP
            } // 170
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 180
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = TEMP
               } // 180
            }
         } // 190

         GO TO 350

         // Use Francis double-shift

         // Note: the Francis double-shift should work with real shifts,
               // but only if the block is at least 3x3.
               // This code may break if this point is reached with
               // a 2x2 block with real eigenvalues.

         } // 200
         if ( IFIRST+1.EQ.ILAST ) {

            // Special case -- 2x2 block with complex eigenvectors

            // Step 1: Standardize, that is, rotate so that

                        // ( B11  0  )
                    // B = (         )  with B11 non-negative.
                        // (  0  B22 )

            dlasv2(T( ILAST-1, ILAST-1 ), T( ILAST-1, ILAST ), T( ILAST, ILAST ), B22, B11, SR, CR, SL, CL );

            if ( B11.LT.ZERO ) {
               CR = -CR
               SR = -SR
               B11 = -B11
               B22 = -B22
            }

            drot(ILASTM+1-IFIRST, H( ILAST-1, ILAST-1 ), LDH, H( ILAST, ILAST-1 ), LDH, CL, SL )             CALL DROT( ILAST+1-IFRSTM, H( IFRSTM, ILAST-1 ), 1, H( IFRSTM, ILAST ), 1, CR, SR );

            IF( ILAST.LT.ILASTM ) CALL DROT( ILASTM-ILAST, T( ILAST-1, ILAST+1 ), LDT, T( ILAST, ILAST+1 ), LDT, CL, SL )             IF( IFRSTM.LT.ILAST-1 ) CALL DROT( IFIRST-IFRSTM, T( IFRSTM, ILAST-1 ), 1, T( IFRSTM, ILAST ), 1, CR, SR )

            IF( ILQ ) CALL DROT( N, Q( 1, ILAST-1 ), 1, Q( 1, ILAST ), 1, CL, SL )             IF( ILZ ) CALL DROT( N, Z( 1, ILAST-1 ), 1, Z( 1, ILAST ), 1, CR, SR )

            T( ILAST-1, ILAST-1 ) = B11
            T( ILAST-1, ILAST ) = ZERO
            T( ILAST, ILAST-1 ) = ZERO
            T( ILAST, ILAST ) = B22

            // If B22 is negative, negate column ILAST

            if ( B22.LT.ZERO ) {
               for (J = IFRSTM; J <= ILAST; J++) { // 210
                  H( J, ILAST ) = -H( J, ILAST )
                  T( J, ILAST ) = -T( J, ILAST )
               } // 210

               if ( ILZ ) {
                  for (J = 1; J <= N; J++) { // 220
                     Z( J, ILAST ) = -Z( J, ILAST )
                  } // 220
               }
               B22 = -B22
            }

            // Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)

            // Recompute shift

            dlag2(H( ILAST-1, ILAST-1 ), LDH, T( ILAST-1, ILAST-1 ), LDT, SAFMIN*SAFETY, S1, TEMP, WR, TEMP2, WI );

            // If standardization has perturbed the shift onto real line,
            // do another (real single-shift) QR step.

            IF( WI.EQ.ZERO ) GO TO 350
            S1INV = ONE / S1

            // Do EISPACK (QZVAL) computation of alpha and beta

            A11 = H( ILAST-1, ILAST-1 )
            A21 = H( ILAST, ILAST-1 )
            A12 = H( ILAST-1, ILAST )
            A22 = H( ILAST, ILAST )

            // Compute complex Givens rotation on right
            // (Assume some element of C = (sA - wB) > unfl )
                             // __
            // (sA - wB) ( CZ   -SZ )
                      // ( SZ    CZ )

            C11R = S1*A11 - WR*B11
            C11I = -WI*B11
            C12 = S1*A12
            C21 = S1*A21
            C22R = S1*A22 - WR*B22
            C22I = -WI*B22

            if ( ABS( C11R )+ABS( C11I )+ABS( C12 ).GT.ABS( C21 )+ ABS( C22R )+ABS( C22I ) ) {
               T1 = DLAPY3( C12, C11R, C11I )
               CZ = C12 / T1
               SZR = -C11R / T1
               SZI = -C11I / T1
            } else {
               CZ = DLAPY2( C22R, C22I )
               if ( CZ.LE.SAFMIN ) {
                  CZ = ZERO
                  SZR = ONE
                  SZI = ZERO
               } else {
                  TEMPR = C22R / CZ
                  TEMPI = C22I / CZ
                  T1 = DLAPY2( CZ, C21 )
                  CZ = CZ / T1
                  SZR = -C21*TEMPR / T1
                  SZI = C21*TEMPI / T1
               }
            }

            // Compute Givens rotation on left

            // (  CQ   SQ )
            // (  __      )  A or B
            // ( -SQ   CQ )

            AN = ABS( A11 ) + ABS( A12 ) + ABS( A21 ) + ABS( A22 )
            BN = ABS( B11 ) + ABS( B22 )
            WABS = ABS( WR ) + ABS( WI )
            if ( S1*AN.GT.WABS*BN ) {
               CQ = CZ*B11
               SQR = SZR*B22
               SQI = -SZI*B22
            } else {
               A1R = CZ*A11 + SZR*A12
               A1I = SZI*A12
               A2R = CZ*A21 + SZR*A22
               A2I = SZI*A22
               CQ = DLAPY2( A1R, A1I )
               if ( CQ.LE.SAFMIN ) {
                  CQ = ZERO
                  SQR = ONE
                  SQI = ZERO
               } else {
                  TEMPR = A1R / CQ
                  TEMPI = A1I / CQ
                  SQR = TEMPR*A2R + TEMPI*A2I
                  SQI = TEMPI*A2R - TEMPR*A2I
               }
            }
            T1 = DLAPY3( CQ, SQR, SQI )
            CQ = CQ / T1
            SQR = SQR / T1
            SQI = SQI / T1

            // Compute diagonal elements of QBZ

            TEMPR = SQR*SZR - SQI*SZI
            TEMPI = SQR*SZI + SQI*SZR
            B1R = CQ*CZ*B11 + TEMPR*B22
            B1I = TEMPI*B22
            B1A = DLAPY2( B1R, B1I )
            B2R = CQ*CZ*B22 + TEMPR*B11
            B2I = -TEMPI*B11
            B2A = DLAPY2( B2R, B2I )

            // Normalize so beta > 0, and Im( alpha1 ) > 0

            BETA( ILAST-1 ) = B1A
            BETA( ILAST ) = B2A
            ALPHAR( ILAST-1 ) = ( WR*B1A )*S1INV
            ALPHAI( ILAST-1 ) = ( WI*B1A )*S1INV
            ALPHAR( ILAST ) = ( WR*B2A )*S1INV
            ALPHAI( ILAST ) = -( WI*B2A )*S1INV

            // Step 3: Go to next block -- exit if finished.

            ILAST = IFIRST - 1
            IF( ILAST.LT.ILO ) GO TO 380

            // Reset counters

            IITER = 0
            ESHIFT = ZERO
            if ( .NOT.ILSCHR ) {
               ILASTM = ILAST
               IF( IFRSTM.GT.ILAST ) IFRSTM = ILO
            }
            GO TO 350
         } else {

            // Usual case: 3x3 or larger block, using Francis implicit
                        // double-shift

                                     // 2
            // Eigenvalue equation is  w  - c w + d = 0,

                                          // -1 2        -1
            // so compute 1st column of  (A B  )  - c A B   + d
            // using the formula in QZIT (from EISPACK)

            // We assume that the block is at least 3x3

            AD11 = ( ASCALE*H( ILAST-1, ILAST-1 ) ) / ( BSCALE*T( ILAST-1, ILAST-1 ) )             AD21 = ( ASCALE*H( ILAST, ILAST-1 ) ) / ( BSCALE*T( ILAST-1, ILAST-1 ) )             AD12 = ( ASCALE*H( ILAST-1, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) )             AD22 = ( ASCALE*H( ILAST, ILAST ) ) / ( BSCALE*T( ILAST, ILAST ) )
            U12 = T( ILAST-1, ILAST ) / T( ILAST, ILAST )
            AD11L = ( ASCALE*H( IFIRST, IFIRST ) ) / ( BSCALE*T( IFIRST, IFIRST ) )             AD21L = ( ASCALE*H( IFIRST+1, IFIRST ) ) / ( BSCALE*T( IFIRST, IFIRST ) )             AD12L = ( ASCALE*H( IFIRST, IFIRST+1 ) ) / ( BSCALE*T( IFIRST+1, IFIRST+1 ) )             AD22L = ( ASCALE*H( IFIRST+1, IFIRST+1 ) ) / ( BSCALE*T( IFIRST+1, IFIRST+1 ) )             AD32L = ( ASCALE*H( IFIRST+2, IFIRST+1 ) ) / ( BSCALE*T( IFIRST+1, IFIRST+1 ) )
            U12L = T( IFIRST, IFIRST+1 ) / T( IFIRST+1, IFIRST+1 )

            V( 1 ) = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 + AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L             V( 2 ) = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )- ( AD22-AD11L )+AD21*U12 )*AD21L
            V( 3 ) = AD32L*AD21L

            ISTART = IFIRST

            dlarfg(3, V( 1 ), V( 2 ), 1, TAU );
            V( 1 ) = ONE

            // Sweep

            for (J = ISTART; J <= ILAST - 2; J++) { // 290

               // All but last elements: use 3x3 Householder transforms.

               // Zero (j-1)st column of A

               if ( J.GT.ISTART ) {
                  V( 1 ) = H( J, J-1 )
                  V( 2 ) = H( J+1, J-1 )
                  V( 3 ) = H( J+2, J-1 )

                  dlarfg(3, H( J, J-1 ), V( 2 ), 1, TAU );
                  V( 1 ) = ONE
                  H( J+1, J-1 ) = ZERO
                  H( J+2, J-1 ) = ZERO
               }

               T2 = TAU*V( 2 )
               T3 = TAU*V( 3 )
               for (JC = J; JC <= ILASTM; JC++) { // 230
                  TEMP = H( J, JC )+V( 2 )*H( J+1, JC )+V( 3 )* H( J+2, JC )
                  H( J, JC ) = H( J, JC ) - TEMP*TAU
                  H( J+1, JC ) = H( J+1, JC ) - TEMP*T2
                  H( J+2, JC ) = H( J+2, JC ) - TEMP*T3
                  TEMP2 = T( J, JC )+V( 2 )*T( J+1, JC )+V( 3 )* T( J+2, JC )
                  T( J, JC ) = T( J, JC ) - TEMP2*TAU
                  T( J+1, JC ) = T( J+1, JC ) - TEMP2*T2
                  T( J+2, JC ) = T( J+2, JC ) - TEMP2*T3
               } // 230
               if ( ILQ ) {
                  for (JR = 1; JR <= N; JR++) { // 240
                     TEMP = Q( JR, J )+V( 2 )*Q( JR, J+1 )+V( 3 )* Q( JR, J+2 )
                     Q( JR, J ) = Q( JR, J ) - TEMP*TAU
                     Q( JR, J+1 ) = Q( JR, J+1 ) - TEMP*T2
                     Q( JR, J+2 ) = Q( JR, J+2 ) - TEMP*T3
                  } // 240
               }

               // Zero j-th column of B (see DLAGBC for details)

               // Swap rows to pivot

               ILPIVT = .FALSE.
               TEMP = MAX( ABS( T( J+1, J+1 ) ), ABS( T( J+1, J+2 ) ) )
               TEMP2 = MAX( ABS( T( J+2, J+1 ) ), ABS( T( J+2, J+2 ) ) )
               if ( MAX( TEMP, TEMP2 ).LT.SAFMIN ) {
                  SCALE = ZERO
                  U1 = ONE
                  U2 = ZERO
                  GO TO 250
               } else if ( TEMP.GE.TEMP2 ) {
                  W11 = T( J+1, J+1 )
                  W21 = T( J+2, J+1 )
                  W12 = T( J+1, J+2 )
                  W22 = T( J+2, J+2 )
                  U1 = T( J+1, J )
                  U2 = T( J+2, J )
               } else {
                  W21 = T( J+1, J+1 )
                  W11 = T( J+2, J+1 )
                  W22 = T( J+1, J+2 )
                  W12 = T( J+2, J+2 )
                  U2 = T( J+1, J )
                  U1 = T( J+2, J )
               }

               // Swap columns if nec.

               if ( ABS( W12 ).GT.ABS( W11 ) ) {
                  ILPIVT = .TRUE.
                  TEMP = W12
                  TEMP2 = W22
                  W12 = W11
                  W22 = W21
                  W11 = TEMP
                  W21 = TEMP2
               }

               // LU-factor

               TEMP = W21 / W11
               U2 = U2 - TEMP*U1
               W22 = W22 - TEMP*W12
               W21 = ZERO

               // Compute SCALE

               SCALE = ONE
               if ( ABS( W22 ).LT.SAFMIN ) {
                  SCALE = ZERO
                  U2 = ONE
                  U1 = -W12 / W11
                  GO TO 250
               }
               IF( ABS( W22 ).LT.ABS( U2 ) ) SCALE = ABS( W22 / U2 )                IF( ABS( W11 ).LT.ABS( U1 ) ) SCALE = MIN( SCALE, ABS( W11 / U1 ) )

               // Solve

               U2 = ( SCALE*U2 ) / W22
               U1 = ( SCALE*U1-W12*U2 ) / W11

               } // 250
               if ( ILPIVT ) {
                  TEMP = U2
                  U2 = U1
                  U1 = TEMP
               }

               // Compute Householder Vector

               T1 = SQRT( SCALE**2+U1**2+U2**2 )
               TAU = ONE + SCALE / T1
               VS = -ONE / ( SCALE+T1 )
               V( 1 ) = ONE
               V( 2 ) = VS*U1
               V( 3 ) = VS*U2

               // Apply transformations from the right.

               T2 = TAU*V(2)
               T3 = TAU*V(3)
               DO 260 JR = IFRSTM, MIN( J+3, ILAST )
                  TEMP = H( JR, J )+V( 2 )*H( JR, J+1 )+V( 3 )* H( JR, J+2 )
                  H( JR, J ) = H( JR, J ) - TEMP*TAU
                  H( JR, J+1 ) = H( JR, J+1 ) - TEMP*T2
                  H( JR, J+2 ) = H( JR, J+2 ) - TEMP*T3
               } // 260
               for (JR = IFRSTM; JR <= J + 2; JR++) { // 270
                  TEMP = T( JR, J )+V( 2 )*T( JR, J+1 )+V( 3 )* T( JR, J+2 )
                  T( JR, J ) = T( JR, J ) - TEMP*TAU
                  T( JR, J+1 ) = T( JR, J+1 ) - TEMP*T2
                  T( JR, J+2 ) = T( JR, J+2 ) - TEMP*T3
               } // 270
               if ( ILZ ) {
                  for (JR = 1; JR <= N; JR++) { // 280
                     TEMP = Z( JR, J )+V( 2 )*Z( JR, J+1 )+V( 3 )* Z( JR, J+2 )
                     Z( JR, J ) = Z( JR, J ) - TEMP*TAU
                     Z( JR, J+1 ) = Z( JR, J+1 ) - TEMP*T2
                     Z( JR, J+2 ) = Z( JR, J+2 ) - TEMP*T3
                  } // 280
               }
               T( J+1, J ) = ZERO
               T( J+2, J ) = ZERO
            } // 290

            // Last elements: Use Givens rotations

            // Rotations from the left

            J = ILAST - 1
            TEMP = H( J, J-1 )
            dlartg(TEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) );
            H( J+1, J-1 ) = ZERO

            for (JC = J; JC <= ILASTM; JC++) { // 300
               TEMP = C*H( J, JC ) + S*H( J+1, JC )
               H( J+1, JC ) = -S*H( J, JC ) + C*H( J+1, JC )
               H( J, JC ) = TEMP
               TEMP2 = C*T( J, JC ) + S*T( J+1, JC )
               T( J+1, JC ) = -S*T( J, JC ) + C*T( J+1, JC )
               T( J, JC ) = TEMP2
            } // 300
            if ( ILQ ) {
               for (JR = 1; JR <= N; JR++) { // 310
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = TEMP
               } // 310
            }

            // Rotations from the right.

            TEMP = T( J+1, J+1 )
            dlartg(TEMP, T( J+1, J ), C, S, T( J+1, J+1 ) );
            T( J+1, J ) = ZERO

            for (JR = IFRSTM; JR <= ILAST; JR++) { // 320
               TEMP = C*H( JR, J+1 ) + S*H( JR, J )
               H( JR, J ) = -S*H( JR, J+1 ) + C*H( JR, J )
               H( JR, J+1 ) = TEMP
            } // 320
            for (JR = IFRSTM; JR <= ILAST - 1; JR++) { // 330
               TEMP = C*T( JR, J+1 ) + S*T( JR, J )
               T( JR, J ) = -S*T( JR, J+1 ) + C*T( JR, J )
               T( JR, J+1 ) = TEMP
            } // 330
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 340
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = TEMP
               } // 340
            }

            // End of Double-Shift code

         }

         GO TO 350

         // End of iteration loop

         } // 350
      } // 360

      // Drop-through = non-convergence

      INFO = ILAST
      GO TO 420

      // Successful completion of all QZ steps

      } // 380

      // Set Eigenvalues 1:ILO-1

      for (J = 1; J <= ILO - 1; J++) { // 410
         if ( T( J, J ).LT.ZERO ) {
            if ( ILSCHR ) {
               for (JR = 1; JR <= J; JR++) { // 390
                  H( JR, J ) = -H( JR, J )
                  T( JR, J ) = -T( JR, J )
               } // 390
            } else {
               H( J, J ) = -H( J, J )
               T( J, J ) = -T( J, J )
            }
            if ( ILZ ) {
               for (JR = 1; JR <= N; JR++) { // 400
                  Z( JR, J ) = -Z( JR, J )
               } // 400
            }
         }
         ALPHAR( J ) = H( J, J )
         ALPHAI( J ) = ZERO
         BETA( J ) = T( J, J )
      } // 410

      // Normal Termination

      INFO = 0

      // Exit (other than argument error) -- return optimal workspace size

      } // 420
      WORK( 1 ) = DBLE( N )
      RETURN

      // End of DHGEQZ

      }

      SUBROUTINE DTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, LWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, LWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      int                NBMIN, NBMAX;
      const              NBMIN = 8, NBMAX = 128 ;
      // ..
      // .. Local Scalars ..
      bool               ALLV, BOTHV, LEFTV, LQUERY, OVER, PAIR, RIGHTV, SOMEV;
      int                I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, IV, MAXWRK, NB, KI2;
      double             BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX, ILAENV;
      double             DDOT, DLAMCH;
      // EXTERNAL LSAME, IDAMAX, ILAENV, DDOT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DLALN2, DSCAL, XERBLA, DGEMM, DLASET, DLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Local Arrays ..
      double             X( 2, 2 );
      int                ISCOMPLEX( NBMAX );
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      BOTHV  = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) || BOTHV
      LEFTV  = LSAME( SIDE, 'L' ) || BOTHV

      ALLV  = LSAME( HOWMNY, 'A' )
      OVER  = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )

      INFO = 0
      NB = ILAENV( 1, 'DTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
      MAXWRK = MAX( 1, N + 2*N*NB )
      WORK(1) = MAXWRK
      LQUERY = ( LWORK == -1 )
      if ( .NOT.RIGHTV && .NOT.LEFTV ) {
         INFO = -1
      } else if ( .NOT.ALLV && .NOT.OVER && .NOT.SOMEV ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDT < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVL < 1 || ( LEFTV && LDVL < N ) ) {
         INFO = -8
      } else if ( LDVR < 1 || ( RIGHTV && LDVR < N ) ) {
         INFO = -10
      } else if ( LWORK < MAX( 1, 3*N ) && .NOT.LQUERY ) {
         INFO = -14
      } else {

         // Set M to the number of columns required to store the selected
         // eigenvectors, standardize the array SELECT if necessary, and
         // test MM.

         if ( SOMEV ) {
            M = 0
            PAIR = false;
            for (J = 1; J <= N; J++) { // 10
               if ( PAIR ) {
                  PAIR = false;
                  SELECT( J ) = false;
               } else {
                  if ( J < N ) {
                     if ( T( J+1, J ) == ZERO ) {
                        IF( SELECT( J ) ) M = M + 1
                     } else {
                        PAIR = true;
                        if ( SELECT( J ) || SELECT( J+1 ) ) {
                           SELECT( J ) = true;
                           M = M + 2
                        }
                     }
                  } else {
                     IF( SELECT( N ) ) M = M + 1
                  }
               }
            } // 10
         } else {
            M = N
         }

         if ( MM < M ) {
            INFO = -11
         }
      }
      if ( INFO != 0 ) {
         xerbla('DTREVC3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      // Use blocked version of back-transformation if sufficient workspace.
      // Zero-out the workspace to avoid potential NaN propagation.

      if ( OVER && LWORK .GE. N + 2*N*NBMIN ) {
         NB = (LWORK - N) / (2*N)
         NB = MIN( NB, NBMAX )
         dlaset('F', N, 1+2*NB, ZERO, ZERO, WORK, N );
      } else {
         NB = 1
      }

      // Set the constants to control overflow.

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      WORK( 1 ) = ZERO
      for (J = 2; J <= N; J++) { // 30
         WORK( J ) = ZERO
         for (I = 1; I <= J - 1; I++) { // 20
            WORK( J ) = WORK( J ) + ABS( T( I, J ) )
         } // 20
      } // 30

      // Index IP is used to specify the real or complex eigenvalue:
        // IP = 0, real eigenvalue,
             // 1, first  of conjugate complex pair: (wr,wi)
            // -1, second of conjugate complex pair: (wr,wi)
        // ISCOMPLEX array stores IP for each column in current block.

      if ( RIGHTV ) {

         // ============================================================
         // Compute right eigenvectors.

         // IV is index of column in current block.
         // For complex right vector, uses IV-1 for real part and IV for complex part.
         // Non-blocked version always uses IV=2;
         // blocked     version starts with IV=NB, goes down to 1 or 2.
         // (Note the "0-th" column is used for 1-norms computed above.)
         IV = 2
         if ( NB.GT.2 ) {
            IV = NB
         }

         IP = 0
         IS = M
         DO 140 KI = N, 1, -1
            if ( IP == -1 ) {
               // previous iteration (ki+1) was second of conjugate pair,
               // so this ki is first of conjugate pair; skip to end of loop
               IP = 1
               GO TO 140
            } else if ( KI == 1 ) {
               // last column, so this ki must be real eigenvalue
               IP = 0
            } else if ( T( KI, KI-1 ) == ZERO ) {
               // zero on sub-diagonal, so this ki is real eigenvalue
               IP = 0
            } else {
               // non-zero on sub-diagonal, so this ki is second of conjugate pair
               IP = -1
            }

            if ( SOMEV ) {
               if ( IP == 0 ) {
                  IF( .NOT.SELECT( KI ) ) GO TO 140
               } else {
                  IF( .NOT.SELECT( KI-1 ) ) GO TO 140
               }
            }

            // Compute the KI-th eigenvalue (WR,WI).

            WR = T( KI, KI )
            WI = ZERO
            if (IP != 0) WI = SQRT( ABS( T( KI, KI-1 ) ) )* SQRT( ABS( T( KI-1, KI ) ) );
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )

            if ( IP == 0 ) {

               // --------------------------------------------------------
               // Real right eigenvector

               WORK( KI + IV*N ) = ONE

               // Form right-hand side.

               for (K = 1; K <= KI - 1; K++) { // 50
                  WORK( K + IV*N ) = -T( K, KI )
               } // 50

               // Solve upper quasi-triangular system:
               // [ T(1:KI-1,1:KI-1) - WR ]*X = SCALE*WORK.

               JNXT = KI - 1
               DO 60 J = KI - 1, 1, -1
                  if (J.GT.JNXT) GO TO 60;
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  if ( J.GT.1 ) {
                     if ( T( J, J-1 ) != ZERO ) {
                        J1   = J - 1
                        JNXT = J - 2
                     }
                  }

                  if ( J1 == J2 ) {

                     // 1-by-1 diagonal block

                     dlaln2( false , 1, 1, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+IV*N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale X(1,1) to avoid overflow when updating
                     // the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        if ( WORK( J ).GT.BIGNUM / XNORM ) {
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        }
                     }

                     // Scale if necessary

                     if (SCALE != ONE) CALL DSCAL( KI, SCALE, WORK( 1+IV*N ), 1 );
                     WORK( J+IV*N ) = X( 1, 1 )

                     // Update right-hand side

                     daxpy(J-1, -X( 1, 1 ), T( 1, J ), 1, WORK( 1+IV*N ), 1 );

                  } else {

                     // 2-by-2 diagonal block

                     dlaln2( false , 2, 1, SMIN, ONE, T( J-1, J-1 ), LDT, ONE, ONE, WORK( J-1+IV*N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale X(1,1) and X(2,1) to avoid overflow when
                     // updating the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        if ( BETA.GT.BIGNUM / XNORM ) {
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 2, 1 ) = X( 2, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        }
                     }

                     // Scale if necessary

                     if (SCALE != ONE) CALL DSCAL( KI, SCALE, WORK( 1+IV*N ), 1 );
                     WORK( J-1+IV*N ) = X( 1, 1 )
                     WORK( J  +IV*N ) = X( 2, 1 )

                     // Update right-hand side

                     daxpy(J-2, -X( 1, 1 ), T( 1, J-1 ), 1, WORK( 1+IV*N ), 1 );
                     daxpy(J-2, -X( 2, 1 ), T( 1, J ), 1, WORK( 1+IV*N ), 1 );
                  }
               } // 60

               // Copy the vector x or Q*x to VR and normalize.

               if ( .NOT.OVER ) {
                  // ------------------------------
                  // no back-transform: copy x to VR and normalize.
                  dcopy(KI, WORK( 1 + IV*N ), 1, VR( 1, IS ), 1 );

                  II = IDAMAX( KI, VR( 1, IS ), 1 )
                  REMAX = ONE / ABS( VR( II, IS ) )
                  dscal(KI, REMAX, VR( 1, IS ), 1 );

                  for (K = KI + 1; K <= N; K++) { // 70
                     VR( K, IS ) = ZERO
                  } // 70

               } else if ( NB == 1 ) {
                  // ------------------------------
                  // version 1: back-transform each vector with GEMV, Q*x.
                  if (KI.GT.1) CALL DGEMV( 'N', N, KI-1, ONE, VR, LDVR, WORK( 1 + IV*N ), 1, WORK( KI + IV*N ), VR( 1, KI ), 1 );

                  II = IDAMAX( N, VR( 1, KI ), 1 )
                  REMAX = ONE / ABS( VR( II, KI ) )
                  dscal(N, REMAX, VR( 1, KI ), 1 );

               } else {
                  // ------------------------------
                  // version 2: back-transform block of vectors with GEMM
                  // zero out below vector
                  for (K = KI + 1; K <= N; K++) {
                     WORK( K + IV*N ) = ZERO
                  }
                  ISCOMPLEX( IV ) = IP
                  // back-transform and normalization is done below
               }
            } else {

               // --------------------------------------------------------
               // Complex right eigenvector.

               // Initial solve
               // [ ( T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I*WI) ]*X = 0.
               // [ ( T(KI,  KI-1) T(KI,  KI) )               ]

               if ( ABS( T( KI-1, KI ) ).GE.ABS( T( KI, KI-1 ) ) ) {
                  WORK( KI-1 + (IV-1)*N ) = ONE
                  WORK( KI   + (IV  )*N ) = WI / T( KI-1, KI )
               } else {
                  WORK( KI-1 + (IV-1)*N ) = -WI / T( KI, KI-1 )
                  WORK( KI   + (IV  )*N ) = ONE
               }
               WORK( KI   + (IV-1)*N ) = ZERO
               WORK( KI-1 + (IV  )*N ) = ZERO

               // Form right-hand side.

               for (K = 1; K <= KI - 2; K++) { // 80
                  WORK( K+(IV-1)*N ) = -WORK( KI-1+(IV-1)*N )*T(K,KI-1)
                  WORK( K+(IV  )*N ) = -WORK( KI  +(IV  )*N )*T(K,KI  )
               } // 80

               // Solve upper quasi-triangular system:
               // [ T(1:KI-2,1:KI-2) - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2)

               JNXT = KI - 2
               DO 90 J = KI - 2, 1, -1
                  if (J.GT.JNXT) GO TO 90;
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  if ( J.GT.1 ) {
                     if ( T( J, J-1 ) != ZERO ) {
                        J1   = J - 1
                        JNXT = J - 2
                     }
                  }

                  if ( J1 == J2 ) {

                     // 1-by-1 diagonal block

                     dlaln2( false , 1, 2, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+(IV-1)*N ), N, WR, WI, X, 2, SCALE, XNORM, IERR );

                     // Scale X(1,1) and X(1,2) to avoid overflow when
                     // updating the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        if ( WORK( J ).GT.BIGNUM / XNORM ) {
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 1, 2 ) = X( 1, 2 ) / XNORM
                           SCALE = SCALE / XNORM
                        }
                     }

                     // Scale if necessary

                     if ( SCALE != ONE ) {
                        dscal(KI, SCALE, WORK( 1+(IV-1)*N ), 1 );
                        dscal(KI, SCALE, WORK( 1+(IV  )*N ), 1 );
                     }
                     WORK( J+(IV-1)*N ) = X( 1, 1 )
                     WORK( J+(IV  )*N ) = X( 1, 2 )

                     // Update the right-hand side

                     daxpy(J-1, -X( 1, 1 ), T( 1, J ), 1, WORK( 1+(IV-1)*N ), 1 );
                     daxpy(J-1, -X( 1, 2 ), T( 1, J ), 1, WORK( 1+(IV  )*N ), 1 );

                  } else {

                     // 2-by-2 diagonal block

                     dlaln2( false , 2, 2, SMIN, ONE, T( J-1, J-1 ), LDT, ONE, ONE, WORK( J-1+(IV-1)*N ), N, WR, WI, X, 2, SCALE, XNORM, IERR );

                     // Scale X to avoid overflow when updating
                     // the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        if ( BETA.GT.BIGNUM / XNORM ) {
                           REC = ONE / XNORM
                           X( 1, 1 ) = X( 1, 1 )*REC
                           X( 1, 2 ) = X( 1, 2 )*REC
                           X( 2, 1 ) = X( 2, 1 )*REC
                           X( 2, 2 ) = X( 2, 2 )*REC
                           SCALE = SCALE*REC
                        }
                     }

                     // Scale if necessary

                     if ( SCALE != ONE ) {
                        dscal(KI, SCALE, WORK( 1+(IV-1)*N ), 1 );
                        dscal(KI, SCALE, WORK( 1+(IV  )*N ), 1 );
                     }
                     WORK( J-1+(IV-1)*N ) = X( 1, 1 )
                     WORK( J  +(IV-1)*N ) = X( 2, 1 )
                     WORK( J-1+(IV  )*N ) = X( 1, 2 )
                     WORK( J  +(IV  )*N ) = X( 2, 2 )

                     // Update the right-hand side

                     daxpy(J-2, -X( 1, 1 ), T( 1, J-1 ), 1, WORK( 1+(IV-1)*N   ), 1 );
                     daxpy(J-2, -X( 2, 1 ), T( 1, J ), 1, WORK( 1+(IV-1)*N   ), 1 );
                     daxpy(J-2, -X( 1, 2 ), T( 1, J-1 ), 1, WORK( 1+(IV  )*N ), 1 );
                     daxpy(J-2, -X( 2, 2 ), T( 1, J ), 1, WORK( 1+(IV  )*N ), 1 );
                  }
               } // 90

               // Copy the vector x or Q*x to VR and normalize.

               if ( .NOT.OVER ) {
                  // ------------------------------
                  // no back-transform: copy x to VR and normalize.
                  dcopy(KI, WORK( 1+(IV-1)*N ), 1, VR(1,IS-1), 1 );
                  dcopy(KI, WORK( 1+(IV  )*N ), 1, VR(1,IS  ), 1 );

                  EMAX = ZERO
                  for (K = 1; K <= KI; K++) { // 100
                     EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+ ABS( VR( K, IS   ) ) )
                  } // 100
                  REMAX = ONE / EMAX
                  dscal(KI, REMAX, VR( 1, IS-1 ), 1 );
                  dscal(KI, REMAX, VR( 1, IS   ), 1 );

                  for (K = KI + 1; K <= N; K++) { // 110
                     VR( K, IS-1 ) = ZERO
                     VR( K, IS   ) = ZERO
                  } // 110

               } else if ( NB == 1 ) {
                  // ------------------------------
                  // version 1: back-transform each vector with GEMV, Q*x.
                  if ( KI.GT.2 ) {
                     dgemv('N', N, KI-2, ONE, VR, LDVR, WORK( 1    + (IV-1)*N ), 1, WORK( KI-1 + (IV-1)*N ), VR(1,KI-1), 1);
                     dgemv('N', N, KI-2, ONE, VR, LDVR, WORK( 1  + (IV)*N ), 1, WORK( KI + (IV)*N ), VR( 1, KI ), 1 );
                  } else {
                     dscal(N, WORK(KI-1+(IV-1)*N), VR(1,KI-1), 1);
                     dscal(N, WORK(KI  +(IV  )*N), VR(1,KI  ), 1);
                  }

                  EMAX = ZERO
                  for (K = 1; K <= N; K++) { // 120
                     EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+ ABS( VR( K, KI   ) ) )
                  } // 120
                  REMAX = ONE / EMAX
                  dscal(N, REMAX, VR( 1, KI-1 ), 1 );
                  dscal(N, REMAX, VR( 1, KI   ), 1 );

               } else {
                  // ------------------------------
                  // version 2: back-transform block of vectors with GEMM
                  // zero out below vector
                  for (K = KI + 1; K <= N; K++) {
                     WORK( K + (IV-1)*N ) = ZERO
                     WORK( K + (IV  )*N ) = ZERO
                  }
                  ISCOMPLEX( IV-1 ) = -IP
                  ISCOMPLEX( IV   ) =  IP
                  IV = IV - 1
                  // back-transform and normalization is done below
               }
            }

            if ( NB.GT.1 ) {
               // --------------------------------------------------------
               // Blocked version of back-transform
               // For complex case, KI2 includes both vectors (KI-1 and KI)
               if ( IP == 0 ) {
                  KI2 = KI
               } else {
                  KI2 = KI - 1
               }

               // Columns IV:NB of work are valid vectors.
               // When the number of vectors stored reaches NB-1 or NB,
               // or if this was last vector, do the GEMM
               if ( (IV.LE.2) || (KI2 == 1) ) {
                  dgemm('N', 'N', N, NB-IV+1, KI2+NB-IV, ONE, VR, LDVR, WORK( 1 + (IV)*N    ), N, ZERO, WORK( 1 + (NB+IV)*N ), N );
                  // normalize vectors
                  for (K = IV; K <= NB; K++) {
                     if ( ISCOMPLEX(K) == 0 ) {
                        // real eigenvector
                        II = IDAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                        REMAX = ONE / ABS( WORK( II + (NB+K)*N ) )
                     } else if ( ISCOMPLEX(K) == 1 ) {
                        // first eigenvector of conjugate pair
                        EMAX = ZERO
                        for (II = 1; II <= N; II++) {
                           EMAX = MAX( EMAX, ABS( WORK( II + (NB+K  )*N ) )+ ABS( WORK( II + (NB+K+1)*N ) ) )
                        }
                        REMAX = ONE / EMAX
                     // else if ISCOMPLEX(K) == -1
                        // second eigenvector of conjugate pair
                        // reuse same REMAX as previous K
                     }
                     dscal(N, REMAX, WORK( 1 + (NB+K)*N ), 1 );
                  }
                  dlacpy('F', N, NB-IV+1, WORK( 1 + (NB+IV)*N ), N, VR( 1, KI2 ), LDVR );
                  IV = NB
               } else {
                  IV = IV - 1
               }
            END IF ! blocked back-transform

            IS = IS - 1
            if (IP != 0) IS = IS - 1;
         } // 140
      }

      if ( LEFTV ) {

         // ============================================================
         // Compute left eigenvectors.

         // IV is index of column in current block.
         // For complex left vector, uses IV for real part and IV+1 for complex part.
         // Non-blocked version always uses IV=1;
         // blocked     version starts with IV=1, goes up to NB-1 or NB.
         // (Note the "0-th" column is used for 1-norms computed above.)
         IV = 1
         IP = 0
         IS = 1
         for (KI = 1; KI <= N; KI++) { // 260
            if ( IP == 1 ) {
               // previous iteration (ki-1) was first of conjugate pair,
               // so this ki is second of conjugate pair; skip to end of loop
               IP = -1
               GO TO 260
            } else if ( KI == N ) {
               // last column, so this ki must be real eigenvalue
               IP = 0
            } else if ( T( KI+1, KI ) == ZERO ) {
               // zero on sub-diagonal, so this ki is real eigenvalue
               IP = 0
            } else {
               // non-zero on sub-diagonal, so this ki is first of conjugate pair
               IP = 1
            }

            if ( SOMEV ) {
               IF( .NOT.SELECT( KI ) ) GO TO 260
            }

            // Compute the KI-th eigenvalue (WR,WI).

            WR = T( KI, KI )
            WI = ZERO
            if (IP != 0) WI = SQRT( ABS( T( KI, KI+1 ) ) )* SQRT( ABS( T( KI+1, KI ) ) );
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )

            if ( IP == 0 ) {

               // --------------------------------------------------------
               // Real left eigenvector

               WORK( KI + IV*N ) = ONE

               // Form right-hand side.

               for (K = KI + 1; K <= N; K++) { // 160
                  WORK( K + IV*N ) = -T( KI, K )
               } // 160

               // Solve transposed quasi-triangular system:
               // [ T(KI+1:N,KI+1:N) - WR ]**T * X = SCALE*WORK

               VMAX = ONE
               VCRIT = BIGNUM

               JNXT = KI + 1
               for (J = KI + 1; J <= N; J++) { // 170
                  if (J < JNXT) GO TO 170;
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  if ( J < N ) {
                     if ( T( J+1, J ) != ZERO ) {
                        J2 = J + 1
                        JNXT = J + 2
                     }
                  }

                  if ( J1 == J2 ) {

                     // 1-by-1 diagonal block

                     // Scale if necessary to avoid overflow when forming
                     // the right-hand side.

                     if ( WORK( J ).GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK( KI+IV*N ), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+IV*N ) = WORK( J+IV*N ) - DDOT( J-KI-1, T( KI+1, J ), 1, WORK( KI+1+IV*N ), 1 )

                     // Solve [ T(J,J) - WR ]**T * X = WORK

                     dlaln2( false , 1, 1, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+IV*N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     if (SCALE != ONE) CALL DSCAL( N-KI+1, SCALE, WORK( KI+IV*N ), 1 );
                     WORK( J+IV*N ) = X( 1, 1 )
                     VMAX = MAX( ABS( WORK( J+IV*N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  } else {

                     // 2-by-2 diagonal block

                     // Scale if necessary to avoid overflow when forming
                     // the right-hand side.

                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     if ( BETA.GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK( KI+IV*N ), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+IV*N ) = WORK( J+IV*N ) - DDOT( J-KI-1, T( KI+1, J ), 1, WORK( KI+1+IV*N ), 1 )

                     WORK( J+1+IV*N ) = WORK( J+1+IV*N ) - DDOT( J-KI-1, T( KI+1, J+1 ), 1, WORK( KI+1+IV*N ), 1 )

                     // Solve
                     // [ T(J,J)-WR   T(J,J+1)      ]**T * X = SCALE*( WORK1 )
                     // [ T(J+1,J)    T(J+1,J+1)-WR ]                ( WORK2 )

                     dlaln2( true , 2, 1, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+IV*N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     if (SCALE != ONE) CALL DSCAL( N-KI+1, SCALE, WORK( KI+IV*N ), 1 );
                     WORK( J  +IV*N ) = X( 1, 1 )
                     WORK( J+1+IV*N ) = X( 2, 1 )

                     VMAX = MAX( ABS( WORK( J  +IV*N ) ), ABS( WORK( J+1+IV*N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  }
               } // 170

               // Copy the vector x or Q*x to VL and normalize.

               if ( .NOT.OVER ) {
                  // ------------------------------
                  // no back-transform: copy x to VL and normalize.
                  dcopy(N-KI+1, WORK( KI + IV*N ), 1, VL( KI, IS ), 1 );

                  II = IDAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                  REMAX = ONE / ABS( VL( II, IS ) )
                  dscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

                  for (K = 1; K <= KI - 1; K++) { // 180
                     VL( K, IS ) = ZERO
                  } // 180

               } else if ( NB == 1 ) {
                  // ------------------------------
                  // version 1: back-transform each vector with GEMV, Q*x.
                  if (KI < N) CALL DGEMV( 'N', N, N-KI, ONE, VL( 1, KI+1 ), LDVL, WORK( KI+1 + IV*N ), 1, WORK( KI   + IV*N ), VL( 1, KI ), 1 );

                  II = IDAMAX( N, VL( 1, KI ), 1 )
                  REMAX = ONE / ABS( VL( II, KI ) )
                  dscal(N, REMAX, VL( 1, KI ), 1 );

               } else {
                  // ------------------------------
                  // version 2: back-transform block of vectors with GEMM
                  // zero out above vector
                  // could go from KI-NV+1 to KI-1
                  for (K = 1; K <= KI - 1; K++) {
                     WORK( K + IV*N ) = ZERO
                  }
                  ISCOMPLEX( IV ) = IP
                  // back-transform and normalization is done below
               }
            } else {

               // --------------------------------------------------------
               // Complex left eigenvector.

               // Initial solve:
               // [ ( T(KI,KI)    T(KI,KI+1)  )**T - (WR - I* WI) ]*X = 0.
               // [ ( T(KI+1,KI) T(KI+1,KI+1) )                   ]

               if ( ABS( T( KI, KI+1 ) ).GE.ABS( T( KI+1, KI ) ) ) {
                  WORK( KI   + (IV  )*N ) = WI / T( KI, KI+1 )
                  WORK( KI+1 + (IV+1)*N ) = ONE
               } else {
                  WORK( KI   + (IV  )*N ) = ONE
                  WORK( KI+1 + (IV+1)*N ) = -WI / T( KI+1, KI )
               }
               WORK( KI+1 + (IV  )*N ) = ZERO
               WORK( KI   + (IV+1)*N ) = ZERO

               // Form right-hand side.

               for (K = KI + 2; K <= N; K++) { // 190
                  WORK( K+(IV  )*N ) = -WORK( KI  +(IV  )*N )*T(KI,  K)
                  WORK( K+(IV+1)*N ) = -WORK( KI+1+(IV+1)*N )*T(KI+1,K)
               } // 190

               // Solve transposed quasi-triangular system:
               // [ T(KI+2:N,KI+2:N)**T - (WR-i*WI) ]*X = WORK1+i*WORK2

               VMAX = ONE
               VCRIT = BIGNUM

               JNXT = KI + 2
               for (J = KI + 2; J <= N; J++) { // 200
                  if (J < JNXT) GO TO 200;
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  if ( J < N ) {
                     if ( T( J+1, J ) != ZERO ) {
                        J2 = J + 1
                        JNXT = J + 2
                     }
                  }

                  if ( J1 == J2 ) {

                     // 1-by-1 diagonal block

                     // Scale if necessary to avoid overflow when
                     // forming the right-hand side elements.

                     if ( WORK( J ).GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK(KI+(IV  )*N), 1 );
                        dscal(N-KI+1, REC, WORK(KI+(IV+1)*N), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+(IV  )*N ) = WORK( J+(IV)*N ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+(IV)*N ), 1 )                      WORK( J+(IV+1)*N ) = WORK( J+(IV+1)*N ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+(IV+1)*N ), 1 )

                     // Solve [ T(J,J)-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2

                     dlaln2( false , 1, 2, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+IV*N ), N, WR, -WI, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     if ( SCALE != ONE ) {
                        dscal(N-KI+1, SCALE, WORK(KI+(IV  )*N), 1);
                        dscal(N-KI+1, SCALE, WORK(KI+(IV+1)*N), 1);
                     }
                     WORK( J+(IV  )*N ) = X( 1, 1 )
                     WORK( J+(IV+1)*N ) = X( 1, 2 )
                     VMAX = MAX( ABS( WORK( J+(IV  )*N ) ), ABS( WORK( J+(IV+1)*N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  } else {

                     // 2-by-2 diagonal block

                     // Scale if necessary to avoid overflow when forming
                     // the right-hand side elements.

                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     if ( BETA.GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK(KI+(IV  )*N), 1 );
                        dscal(N-KI+1, REC, WORK(KI+(IV+1)*N), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J  +(IV  )*N ) = WORK( J+(IV)*N ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+(IV)*N ), 1 )

                     WORK( J  +(IV+1)*N ) = WORK( J+(IV+1)*N ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+(IV+1)*N ), 1 )

                     WORK( J+1+(IV  )*N ) = WORK( J+1+(IV)*N ) - DDOT( J-KI-2, T( KI+2, J+1 ), 1, WORK( KI+2+(IV)*N ), 1 )

                     WORK( J+1+(IV+1)*N ) = WORK( J+1+(IV+1)*N ) - DDOT( J-KI-2, T( KI+2, J+1 ), 1, WORK( KI+2+(IV+1)*N ), 1 )

                     // Solve 2-by-2 complex linear equation
                     // [ (T(j,j)   T(j,j+1)  )**T - (wr-i*wi)*I ]*X = SCALE*B
                     // [ (T(j+1,j) T(j+1,j+1))                  ]

                     dlaln2( true , 2, 2, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+IV*N ), N, WR, -WI, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     if ( SCALE != ONE ) {
                        dscal(N-KI+1, SCALE, WORK(KI+(IV  )*N), 1);
                        dscal(N-KI+1, SCALE, WORK(KI+(IV+1)*N), 1);
                     }
                     WORK( J  +(IV  )*N ) = X( 1, 1 )
                     WORK( J  +(IV+1)*N ) = X( 1, 2 )
                     WORK( J+1+(IV  )*N ) = X( 2, 1 )
                     WORK( J+1+(IV+1)*N ) = X( 2, 2 )
                     VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ), ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  }
               } // 200

               // Copy the vector x or Q*x to VL and normalize.

               if ( .NOT.OVER ) {
                  // ------------------------------
                  // no back-transform: copy x to VL and normalize.
                  dcopy(N-KI+1, WORK( KI + (IV  )*N ), 1, VL( KI, IS   ), 1 );
                  dcopy(N-KI+1, WORK( KI + (IV+1)*N ), 1, VL( KI, IS+1 ), 1 );

                  EMAX = ZERO
                  for (K = KI; K <= N; K++) { // 220
                     EMAX = MAX( EMAX, ABS( VL( K, IS   ) )+ ABS( VL( K, IS+1 ) ) )
                  } // 220
                  REMAX = ONE / EMAX
                  dscal(N-KI+1, REMAX, VL( KI, IS   ), 1 );
                  dscal(N-KI+1, REMAX, VL( KI, IS+1 ), 1 );

                  for (K = 1; K <= KI - 1; K++) { // 230
                     VL( K, IS   ) = ZERO
                     VL( K, IS+1 ) = ZERO
                  } // 230

               } else if ( NB == 1 ) {
                  // ------------------------------
                  // version 1: back-transform each vector with GEMV, Q*x.
                  if ( KI < N-1 ) {
                     dgemv('N', N, N-KI-1, ONE, VL( 1, KI+2 ), LDVL, WORK( KI+2 + (IV)*N ), 1, WORK( KI   + (IV)*N ), VL( 1, KI ), 1 );
                     dgemv('N', N, N-KI-1, ONE, VL( 1, KI+2 ), LDVL, WORK( KI+2 + (IV+1)*N ), 1, WORK( KI+1 + (IV+1)*N ), VL( 1, KI+1 ), 1 );
                  } else {
                     dscal(N, WORK(KI+  (IV  )*N), VL(1, KI  ), 1);
                     dscal(N, WORK(KI+1+(IV+1)*N), VL(1, KI+1), 1);
                  }

                  EMAX = ZERO
                  for (K = 1; K <= N; K++) { // 240
                     EMAX = MAX( EMAX, ABS( VL( K, KI   ) )+ ABS( VL( K, KI+1 ) ) )
                  } // 240
                  REMAX = ONE / EMAX
                  dscal(N, REMAX, VL( 1, KI   ), 1 );
                  dscal(N, REMAX, VL( 1, KI+1 ), 1 );

               } else {
                  // ------------------------------
                  // version 2: back-transform block of vectors with GEMM
                  // zero out above vector
                  // could go from KI-NV+1 to KI-1
                  for (K = 1; K <= KI - 1; K++) {
                     WORK( K + (IV  )*N ) = ZERO
                     WORK( K + (IV+1)*N ) = ZERO
                  }
                  ISCOMPLEX( IV   ) =  IP
                  ISCOMPLEX( IV+1 ) = -IP
                  IV = IV + 1
                  // back-transform and normalization is done below
               }
            }

            if ( NB.GT.1 ) {
               // --------------------------------------------------------
               // Blocked version of back-transform
               // For complex case, KI2 includes both vectors (KI and KI+1)
               if ( IP == 0 ) {
                  KI2 = KI
               } else {
                  KI2 = KI + 1
               }

               // Columns 1:IV of work are valid vectors.
               // When the number of vectors stored reaches NB-1 or NB,
               // or if this was last vector, do the GEMM
               if ( (IV.GE.NB-1) || (KI2 == N) ) {
                  dgemm('N', 'N', N, IV, N-KI2+IV, ONE, VL( 1, KI2-IV+1 ), LDVL, WORK( KI2-IV+1 + (1)*N ), N, ZERO, WORK( 1 + (NB+1)*N ), N );
                  // normalize vectors
                  for (K = 1; K <= IV; K++) {
                     if ( ISCOMPLEX(K) == 0) {
                        // real eigenvector
                        II = IDAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                        REMAX = ONE / ABS( WORK( II + (NB+K)*N ) )
                     } else if ( ISCOMPLEX(K) == 1) {
                        // first eigenvector of conjugate pair
                        EMAX = ZERO
                        for (II = 1; II <= N; II++) {
                           EMAX = MAX( EMAX, ABS( WORK( II + (NB+K  )*N ) )+ ABS( WORK( II + (NB+K+1)*N ) ) )
                        }
                        REMAX = ONE / EMAX
                     // else if ISCOMPLEX(K) == -1
                        // second eigenvector of conjugate pair
                        // reuse same REMAX as previous K
                     }
                     dscal(N, REMAX, WORK( 1 + (NB+K)*N ), 1 );
                  }
                  dlacpy('F', N, IV, WORK( 1 + (NB+1)*N ), N, VL( 1, KI2-IV+1 ), LDVL );
                  IV = 1
               } else {
                  IV = IV + 1
               }
            END IF ! blocked back-transform

            IS = IS + 1
            if (IP != 0) IS = IS + 1;
         } // 260
      }

      RETURN

      // End of DTREVC3

      }

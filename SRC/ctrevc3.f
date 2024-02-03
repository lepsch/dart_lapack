      SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      REAL               RWORK( * )
      COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE  = ( 1.0E+0, 0.0E+0 ) ;
      int                NBMIN, NBMAX;
      const              NBMIN = 8, NBMAX = 128 ;
      // ..
      // .. Local Scalars ..
      bool               ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV;
      int                I, II, IS, J, K, KI, IV, MAXWRK, NB;
      REAL               OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX            CDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV, ICAMAX;
      REAL               SLAMCH, SCASUM, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, ICAMAX, SLAMCH, SCASUM, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CCOPY, CLASET, CSSCAL, CGEMM, CGEMV, CLATRS, CLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, CMPLX, CONJG, AIMAG, MAX
      // ..
      // .. Statement Functions ..
      REAL   CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      BOTHV  = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV  = LSAME( SIDE, 'L' ) .OR. BOTHV

      ALLV  = LSAME( HOWMNY, 'A' )
      OVER  = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )

      // Set M to the number of columns required to store the selected
      // eigenvectors.

      if ( SOMEV ) {
         M = 0
         for (J = 1; J <= N; J++) { // 10
            IF( SELECT( J ) ) M = M + 1
         } // 10
      } else {
         M = N
      }

      INFO = 0
      NB = ILAENV( 1, 'CTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
      MAXWRK = MAX( 1, N + 2*N*NB )
      WORK(1) = SROUNDUP_LWORK(MAXWRK)
      RWORK(1) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
      if ( .NOT.RIGHTV .AND. .NOT.LEFTV ) {
         INFO = -1
      } else if ( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) {
         INFO = -8
      } else if ( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) {
         INFO = -10
      } else if ( MM.LT.M ) {
         INFO = -11
      } else if ( LWORK.LT.MAX( 1, 2*N ) .AND. .NOT.LQUERY ) {
         INFO = -14
      } else if ( LRWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
         INFO = -16
      }
      if ( INFO.NE.0 ) {
         xerbla('CTREVC3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible.

      IF( N.EQ.0 ) RETURN

      // Use blocked version of back-transformation if sufficient workspace.
      // Zero-out the workspace to avoid potential NaN propagation.

      if ( OVER .AND. LWORK .GE. N + 2*N*NBMIN ) {
         NB = (LWORK - N) / (2*N)
         NB = MIN( NB, NBMAX )
         claset('F', N, 1+2*NB, CZERO, CZERO, WORK, N );
      } else {
         NB = 1
      }

      // Set the constants to control overflow.

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )

      // Store the diagonal elements of T in working array WORK.

      for (I = 1; I <= N; I++) { // 20
         WORK( I ) = T( I, I )
      } // 20

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      RWORK( 1 ) = ZERO
      for (J = 2; J <= N; J++) { // 30
         RWORK( J ) = SCASUM( J-1, T( 1, J ), 1 )
      } // 30

      if ( RIGHTV ) {

         // ============================================================
         // Compute right eigenvectors.

         // IV is index of column in current block.
         // Non-blocked version always uses IV=NB=1;
         // blocked     version starts with IV=NB, goes down to 1.
         // (Note the "0-th" column is used to store the original diagonal.)
         IV = NB
         IS = M
         DO 80 KI = N, 1, -1
            if ( SOMEV ) {
               IF( .NOT.SELECT( KI ) ) GO TO 80
            }
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )

            // --------------------------------------------------------
            // Complex right eigenvector

            WORK( KI + IV*N ) = CONE

            // Form right-hand side.

            DO 40 K = 1, KI - 1
               WORK( K + IV*N ) = -T( K, KI )
            } // 40

            // Solve upper triangular system:
            // [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.

            DO 50 K = 1, KI - 1
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) T( K, K ) = SMIN
            } // 50

            if ( KI.GT.1 ) {
               clatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI-1, T, LDT, WORK( 1 + IV*N ), SCALE, RWORK, INFO );
               WORK( KI + IV*N ) = SCALE
            }

            // Copy the vector x or Q*x to VR and normalize.

            if ( .NOT.OVER ) {
               // ------------------------------
               // no back-transform: copy x to VR and normalize.
               ccopy(KI, WORK( 1 + IV*N ), 1, VR( 1, IS ), 1 );

               II = ICAMAX( KI, VR( 1, IS ), 1 )
               REMAX = ONE / CABS1( VR( II, IS ) )
               csscal(KI, REMAX, VR( 1, IS ), 1 );

               DO 60 K = KI + 1, N
                  VR( K, IS ) = CZERO
               } // 60

            } else if ( NB.EQ.1 ) {
               // ------------------------------
               // version 1: back-transform each vector with GEMV, Q*x.
               IF( KI.GT.1 ) CALL CGEMV( 'N', N, KI-1, CONE, VR, LDVR, WORK( 1 + IV*N ), 1, CMPLX( SCALE ), VR( 1, KI ), 1 )

               II = ICAMAX( N, VR( 1, KI ), 1 )
               REMAX = ONE / CABS1( VR( II, KI ) )
               csscal(N, REMAX, VR( 1, KI ), 1 );

            } else {
               // ------------------------------
               // version 2: back-transform block of vectors with GEMM
               // zero out below vector
               DO K = KI + 1, N
                  WORK( K + IV*N ) = CZERO
               END DO

               // Columns IV:NB of work are valid vectors.
               // When the number of vectors stored reaches NB,
               // or if this was last vector, do the GEMM
               if ( (IV.EQ.1) .OR. (KI.EQ.1) ) {
                  cgemm('N', 'N', N, NB-IV+1, KI+NB-IV, CONE, VR, LDVR, WORK( 1 + (IV)*N    ), N, CZERO, WORK( 1 + (NB+IV)*N ), N );
                  // normalize vectors
                  for (K = IV; K <= NB; K++) {
                     II = ICAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) )
                     csscal(N, REMAX, WORK( 1 + (NB+K)*N ), 1 );
                  END DO
                  clacpy('F', N, NB-IV+1, WORK( 1 + (NB+IV)*N ), N, VR( 1, KI ), LDVR );
                  IV = NB
               } else {
                  IV = IV - 1
               }
            }

            // Restore the original diagonal elements of T.

            DO 70 K = 1, KI - 1
               T( K, K ) = WORK( K )
            } // 70

            IS = IS - 1
         } // 80
      }

      if ( LEFTV ) {

         // ============================================================
         // Compute left eigenvectors.

         // IV is index of column in current block.
         // Non-blocked version always uses IV=1;
         // blocked     version starts with IV=1, goes up to NB.
         // (Note the "0-th" column is used to store the original diagonal.)
         IV = 1
         IS = 1
         for (KI = 1; KI <= N; KI++) { // 130

            if ( SOMEV ) {
               IF( .NOT.SELECT( KI ) ) GO TO 130
            }
            SMIN = MAX( ULP*( CABS1( T( KI, KI ) ) ), SMLNUM )

            // --------------------------------------------------------
            // Complex left eigenvector

            WORK( KI + IV*N ) = CONE

            // Form right-hand side.

            DO 90 K = KI + 1, N
               WORK( K + IV*N ) = -CONJG( T( KI, K ) )
            } // 90

            // Solve conjugate-transposed triangular system:
            // [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.

            DO 100 K = KI + 1, N
               T( K, K ) = T( K, K ) - T( KI, KI )
               IF( CABS1( T( K, K ) ).LT.SMIN ) T( K, K ) = SMIN
            } // 100

            if ( KI.LT.N ) {
               clatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N-KI, T( KI+1, KI+1 ), LDT, WORK( KI+1 + IV*N ), SCALE, RWORK, INFO );
               WORK( KI + IV*N ) = SCALE
            }

            // Copy the vector x or Q*x to VL and normalize.

            if ( .NOT.OVER ) {
               // ------------------------------
               // no back-transform: copy x to VL and normalize.
               ccopy(N-KI+1, WORK( KI + IV*N ), 1, VL(KI,IS), 1 );

               II = ICAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
               REMAX = ONE / CABS1( VL( II, IS ) )
               csscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

               DO 110 K = 1, KI - 1
                  VL( K, IS ) = CZERO
               } // 110

            } else if ( NB.EQ.1 ) {
               // ------------------------------
               // version 1: back-transform each vector with GEMV, Q*x.
               IF( KI.LT.N ) CALL CGEMV( 'N', N, N-KI, CONE, VL( 1, KI+1 ), LDVL, WORK( KI+1 + IV*N ), 1, CMPLX( SCALE ), VL( 1, KI ), 1 )

               II = ICAMAX( N, VL( 1, KI ), 1 )
               REMAX = ONE / CABS1( VL( II, KI ) )
               csscal(N, REMAX, VL( 1, KI ), 1 );

            } else {
               // ------------------------------
               // version 2: back-transform block of vectors with GEMM
               // zero out above vector
               // could go from KI-NV+1 to KI-1
               DO K = 1, KI - 1
                  WORK( K + IV*N ) = CZERO
               END DO

               // Columns 1:IV of work are valid vectors.
               // When the number of vectors stored reaches NB,
               // or if this was last vector, do the GEMM
               if ( (IV.EQ.NB) .OR. (KI.EQ.N) ) {
                  cgemm('N', 'N', N, IV, N-KI+IV, CONE, VL( 1, KI-IV+1 ), LDVL, WORK( KI-IV+1 + (1)*N ), N, CZERO, WORK( 1 + (NB+1)*N ), N );
                  // normalize vectors
                  for (K = 1; K <= IV; K++) {
                     II = ICAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                     REMAX = ONE / CABS1( WORK( II + (NB+K)*N ) )
                     csscal(N, REMAX, WORK( 1 + (NB+K)*N ), 1 );
                  END DO
                  clacpy('F', N, IV, WORK( 1 + (NB+1)*N ), N, VL( 1, KI-IV+1 ), LDVL );
                  IV = 1
               } else {
                  IV = IV + 1
               }
            }

            // Restore the original diagonal elements of T.

            DO 120 K = KI + 1, N
               T( K, K ) = WORK( K )
            } // 120

            IS = IS + 1
         } // 130
      }

      RETURN

      // End of CTREVC3

      }

// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void strsna(final int JOB, final int HOWMNY, final int SELECT, final int N, final Matrix<double> T_, final int LDT, final Matrix<double> VL_, final int LDVL, final Matrix<double> VR_, final int LDVR, final int S, final int SEP, final int MM, final int M, final Matrix<double> WORK_, final int LDWORK, final Array<int> IWORK_, final Box<int> INFO,) {
  final T = T_.dim();
  final VL = VL_.dim();
  final VR = VR_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             HOWMNY, JOB;
      int                INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N;
      bool               SELECT( * );
      int                IWORK( * );
      double               S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( LDWORK, * );
      // ..

      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      bool               PAIR, SOMCON, WANTBH, WANTS, WANTSP;
      int                I, IERR, IFST, ILST, J, K, KASE, KS, N2, NN;
      double               BIGNUM, COND, CS, DELTA, DUMM, EPS, EST, LNRM, MU, PROD, PROD1, PROD2, RNRM, SCALE, SMLNUM, SN;
      int                ISAVE( 3 );
      double               DUMMY( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SDOT, SLAMCH, SLAPY2, SNRM2;
      // EXTERNAL lsame, SDOT, SLAMCH, SLAPY2, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLACPY, SLAQTR, STREXC, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT

      // Decode and test the input parameters

      WANTBH = lsame( JOB, 'B' );
      WANTS = lsame( JOB, 'E' ) || WANTBH;
      WANTSP = lsame( JOB, 'V' ) || WANTBH;

      SOMCON = lsame( HOWMNY, 'S' );

      INFO = 0;
      if ( !WANTS && !WANTSP ) {
         INFO = -1;
      } else if ( !lsame( HOWMNY, 'A' ) && !SOMCON ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDVL < 1 || ( WANTS && LDVL < N ) ) {
         INFO = -8;
      } else if ( LDVR < 1 || ( WANTS && LDVR < N ) ) {
         INFO = -10;
      } else {

         // Set M to the number of eigenpairs for which condition numbers
         // are required, and test MM.

         if ( SOMCON ) {
            M = 0;
            PAIR = false;
            for (K = 1; K <= N; K++) { // 10
               if ( PAIR ) {
                  PAIR = false;
               } else {
                  if ( K < N ) {
                     if ( T( K+1, K ) == ZERO ) {
                        if( SELECT( K ) ) M = M + 1;
                     } else {
                        PAIR = true;
                        if( SELECT( K ) || SELECT( K+1 ) ) M = M + 2;
                     }
                  } else {
                     if( SELECT( N ) ) M = M + 1;
                  }
               }
            } // 10
         } else {
            M = N;
         }

         if ( MM < M ) {
            INFO = -13;
         } else if ( LDWORK < 1 || ( WANTSP && LDWORK < N ) ) {
            INFO = -16;
         }
      }
      if ( INFO != 0 ) {
         xerbla('STRSNA', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( N == 1 ) {
         if ( SOMCON ) {
            if( !SELECT( 1 ) ) return;
         }
         if (WANTS) S( 1 ) = ONE;
         IF[WANTSP ) SEP( 1] = ( T( 1, 1 ) ).abs();
         return;
      }

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      KS = 0;
      PAIR = false;
      for (K = 1; K <= N; K++) { // 60

         // Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.

         if ( PAIR ) {
            PAIR = false;
            GO TO 60;
         } else {
            if (K < N) PAIR = T( K+1, K ) != ZERO;
         }

         // Determine whether condition numbers are required for the k-th
         // eigenpair.

         if ( SOMCON ) {
            if ( PAIR ) {
               if( !SELECT( K ) && !SELECT( K+1 ) ) GO TO 60;
            } else {
               if( !SELECT( K ) ) GO TO 60;
            }
         }

         KS = KS + 1;

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            if ( !PAIR ) {

               // Real eigenvalue.

               PROD = SDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 );
               RNRM = SNRM2( N, VR( 1, KS ), 1 );
               LNRM = SNRM2( N, VL( 1, KS ), 1 );
               S[KS] = ( PROD ).abs() / ( RNRM*LNRM );
            } else {

               // Complex eigenvalue.

               PROD1 = SDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 );
               PROD1 = PROD1 + SDOT( N, VR( 1, KS+1 ), 1, VL( 1, KS+1 ), 1 );
               PROD2 = SDOT( N, VL( 1, KS ), 1, VR( 1, KS+1 ), 1 );
               PROD2 = PROD2 - SDOT( N, VL( 1, KS+1 ), 1, VR( 1, KS ), 1 )                RNRM = SLAPY2( SNRM2( N, VR( 1, KS ), 1 ), SNRM2( N, VR( 1, KS+1 ), 1 ) )                LNRM = SLAPY2( SNRM2( N, VL( 1, KS ), 1 ), SNRM2( N, VL( 1, KS+1 ), 1 ) );
               COND = SLAPY2( PROD1, PROD2 ) / ( RNRM*LNRM );
               S[KS] = COND;
               S[KS+1] = COND;
            }
         }

         if ( WANTSP ) {

            // Estimate the reciprocal condition number of the k-th
            // eigenvector.

            // Copy the matrix T to the array WORK and swap the diagonal
            // block beginning at T(k,k) to the (1,1) position.

            slacpy('Full', N, N, T, LDT, WORK, LDWORK );
            IFST = K;
            ILST = 1;
            strexc('No Q', N, WORK, LDWORK, DUMMY, 1, IFST, ILST, WORK( 1, N+1 ), IERR );

            if ( IERR == 1 || IERR == 2 ) {

               // Could not swap because blocks not well separated

               SCALE = ONE;
               EST = BIGNUM;
            } else {

               // Reordering successful

               if ( WORK( 2, 1 ) == ZERO ) {

                  // Form C = T22 - lambda*I in WORK(2:N,2:N).

                  for (I = 2; I <= N; I++) { // 20
                     WORK[I][I] = WORK( I, I ) - WORK( 1, 1 );
                  } // 20
                  N2 = 1;
                  NN = N - 1;
               } else {

                  // Triangularize the 2 by 2 block by unitary
                  // transformation U = [  cs   i*ss ]
                  //                    [ i*ss   cs  ].
                  // such that the (1,1) position of WORK is complex
                  // eigenvalue lambda with positive imaginary part. (2,2)
                  // position of WORK is the complex eigenvalue lambda
                  // with negative imaginary  part.

                  MU = sqrt( ( WORK( 1, 2 ) ).abs() )* sqrt( ( WORK( 2, 1 ) ).abs() );
                  DELTA = SLAPY2( MU, WORK( 2, 1 ) );
                  CS = MU / DELTA;
                  SN = -WORK( 2, 1 ) / DELTA;

                  // Form

                  // C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
                  //                          [   mu                     ]
                  //                          [         ..               ]
                  //                          [             ..           ]
                  //                          [                  mu      ]
                  // where C**T is transpose of matrix C,
                  // and RWORK is stored starting in the N+1-st column of
                  // WORK.

                  for (J = 3; J <= N; J++) { // 30
                     WORK[2][J] = CS*WORK( 2, J );
                     WORK[J][J] = WORK( J, J ) - WORK( 1, 1 );
                  } // 30
                  WORK[2][2] = ZERO;

                  WORK[1][N+1] = TWO*MU;
                  for (I = 2; I <= N - 1; I++) { // 40
                     WORK[I][N+1] = SN*WORK( 1, I+1 );
                  } // 40
                  N2 = 2;
                  NN = 2*( N-1 );
               }

               // Estimate norm(inv(C**T))

               EST = ZERO;
               KASE = 0;
               } // 50
               slacn2(NN, WORK( 1, N+2 ), WORK( 1, N+4 ), IWORK, EST, KASE, ISAVE );
               if ( KASE != 0 ) {
                  if ( KASE == 1 ) {
                     if ( N2 == 1 ) {

                        // Real eigenvalue: solve C**T*x = scale*c.

                        slaqtr( true , true , N-1, WORK( 2, 2 ), LDWORK, DUMMY, DUMM, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR );
                     } else {

                        // Complex eigenvalue: solve
                        // C**T*(p+iq) = scale*(c+id) in real arithmetic.

                        slaqtr( true , false , N-1, WORK( 2, 2 ), LDWORK, WORK( 1, N+1 ), MU, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR );
                     }
                  } else {
                     if ( N2 == 1 ) {

                        // Real eigenvalue: solve C*x = scale*c.

                        slaqtr( false , true , N-1, WORK( 2, 2 ), LDWORK, DUMMY, DUMM, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR );
                     } else {

                        // Complex eigenvalue: solve
                        // C*(p+iq) = scale*(c+id) in real arithmetic.

                        slaqtr( false , false , N-1, WORK( 2, 2 ), LDWORK, WORK( 1, N+1 ), MU, SCALE, WORK( 1, N+4 ), WORK( 1, N+6 ), IERR );

                     }
                  }

                  GO TO 50;
               }
            }

            SEP[KS] = SCALE / max( EST, SMLNUM );
            if (PAIR) SEP( KS+1 ) = SEP( KS );
         }

         if (PAIR) KS = KS + 1;

      } // 60
      }

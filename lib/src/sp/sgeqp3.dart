// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sgeqp3(final int M, final int N, final Matrix<double> A_, final int LDA, final int JPVT, final int TAU, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LWORK, M, N;
      int                JPVT( * );
      double               A( LDA, * ), TAU( * ), WORK( * );
      // ..

      int                INB, INBMIN, IXOVER;
      const              INB = 1, INBMIN = 2, IXOVER = 3 ;
      bool               LQUERY;
      int                FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB, NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SLAQP2, SLAQPS, SORMQR, SSWAP, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SNRM2, SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SNRM2, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // Test input arguments
// ====================

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }

      if ( INFO == 0 ) {
         MINMN = min( M, N );
         if ( MINMN == 0 ) {
            IWS = 1;
            LWKOPT = 1;
         } else {
            IWS = 3*N + 1;
            NB = ilaenv( INB, 'SGEQRF', ' ', M, N, -1, -1 );
            LWKOPT = 2*N + ( N + 1 )*NB;
         }
         WORK[1] = SROUNDUP_LWORK(LWKOPT);

         if ( ( LWORK < IWS ) && !LQUERY ) {
            INFO = -8;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGEQP3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Move initial columns up front.

      NFXD = 1;
      for (J = 1; J <= N; J++) { // 10
         if ( JPVT( J ) != 0 ) {
            if ( J != NFXD ) {
               sswap(M, A( 1, J ), 1, A( 1, NFXD ), 1 );
               JPVT[J] = JPVT( NFXD );
               JPVT[NFXD] = J;
            } else {
               JPVT[J] = J;
            }
            NFXD = NFXD + 1;
         } else {
            JPVT[J] = J;
         }
      } // 10
      NFXD = NFXD - 1;

      // Factorize fixed columns
// =======================

      // Compute the QR factorization of fixed columns and update
      // remaining columns.

      if ( NFXD > 0 ) {
         NA = min( M, NFXD );
// CC      CALL SGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         sgeqrf(M, NA, A, LDA, TAU, WORK, LWORK, INFO );
         IWS = max( IWS, INT( WORK( 1 ) ) );
         if ( NA < N ) {
// CC         CALL SORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
// CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            sormqr('Left', 'Transpose', M, N-NA, NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, LWORK, INFO );
            IWS = max( IWS, INT( WORK( 1 ) ) );
         }
      }

      // Factorize free columns
// ======================

      if ( NFXD < MINMN ) {

         SM = M - NFXD;
         SN = N - NFXD;
         SMINMN = MINMN - NFXD;

         // Determine the block size.

         NB = ilaenv( INB, 'SGEQRF', ' ', SM, SN, -1, -1 );
         NBMIN = 2;
         NX = 0;

         if ( ( NB > 1 ) && ( NB < SMINMN ) ) {

            // Determine when to cross over from blocked to unblocked code.

            NX = max( 0, ilaenv( IXOVER, 'SGEQRF', ' ', SM, SN, -1, -1 ) );


            if ( NX < SMINMN ) {

               // Determine if workspace is large enough for blocked code.

               MINWS = 2*SN + ( SN+1 )*NB;
               IWS = max( IWS, MINWS );
               if ( LWORK < MINWS ) {

                  // Not enough workspace to use optimal NB: Reduce NB and
                  // determine the minimum value of NB.

                  NB = ( LWORK-2*SN ) / ( SN+1 );
                  NBMIN = max( 2, ilaenv( INBMIN, 'SGEQRF', ' ', SM, SN, -1, -1 ) );


               }
            }
         }

         // Initialize partial column norms. The first N elements of work
         // store the exact column norms.

         for (J = NFXD + 1; J <= N; J++) { // 20
            WORK[J] = SNRM2( SM, A( NFXD+1, J ), 1 );
            WORK[N+J] = WORK( J );
         } // 20

         if ( ( NB >= NBMIN ) && ( NB < SMINMN ) && ( NX < SMINMN ) ) {

            // Use blocked code initially.

            J = NFXD + 1;

            // Compute factorization: while loop.


            TOPBMN = MINMN - NX;
            } // 30
            if ( J <= TOPBMN ) {
               JB = min( NB, TOPBMN-J+1 );

               // Factorize JB columns among columns J:N.

               slaqps(M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA, JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ), WORK( 2*N+1 ), WORK( 2*N+JB+1 ), N-J+1 );

               J = J + FJB;
               GO TO 30;
            }
         } else {
            J = NFXD + 1;
         }

         // Use unblocked code to factor the last or only block.


         if (J <= MINMN) slaqp2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ), WORK( 2*N+1 ) );

      }

      WORK[1] = SROUNDUP_LWORK(IWS);
      }

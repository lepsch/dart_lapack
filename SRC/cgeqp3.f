      SUBROUTINE CGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               RWORK( * );
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                INB, INBMIN, IXOVER;
      const              INB = 1, INBMIN = 2, IXOVER = 3 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB, NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CLAQP2, CLAQPS, CSWAP, CUNMQR, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SCNRM2;
      // EXTERNAL ILAENV, SCNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

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
            IWS = N + 1;
            NB = ILAENV( INB, 'CGEQRF', ' ', M, N, -1, -1 );
            LWKOPT = ( N + 1 )*NB;
         }
         WORK( 1 ) = CMPLX( LWKOPT );

         if ( ( LWORK < IWS ) && !LQUERY ) {
            INFO = -8;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGEQP3', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Move initial columns up front.

      NFXD = 1;
      for (J = 1; J <= N; J++) { // 10
         if ( JPVT( J ) != 0 ) {
            if ( J != NFXD ) {
               cswap(M, A( 1, J ), 1, A( 1, NFXD ), 1 );
               JPVT( J ) = JPVT( NFXD );
               JPVT( NFXD ) = J;
            } else {
               JPVT( J ) = J;
            }
            NFXD = NFXD + 1;
         } else {
            JPVT( J ) = J;
         }
      } // 10
      NFXD = NFXD - 1;

      // Factorize fixed columns
// =======================

      // Compute the QR factorization of fixed columns and update
      // remaining columns.

      if ( NFXD > 0 ) {
         NA = min( M, NFXD );
// CC      CALL CGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         cgeqrf(M, NA, A, LDA, TAU, WORK, LWORK, INFO );
         IWS = max( IWS, INT( WORK( 1 ) ) );
         if ( NA < N ) {
// CC         CALL CUNM2R( 'Left', 'Conjugate Transpose', M, N-NA,
// CC  $                   NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK,
// CC  $                   INFO )
            cunmqr('Left', 'Conjugate Transpose', M, N-NA, NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, LWORK, INFO );
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

         NB = ILAENV( INB, 'CGEQRF', ' ', SM, SN, -1, -1 );
         NBMIN = 2;
         NX = 0;

         if ( ( NB > 1 ) && ( NB < SMINMN ) ) {

            // Determine when to cross over from blocked to unblocked code.

            NX = max( 0, ILAENV( IXOVER, 'CGEQRF', ' ', SM, SN, -1, -1 ) );


            if ( NX < SMINMN ) {

               // Determine if workspace is large enough for blocked code.

               MINWS = ( SN+1 )*NB;
               IWS = max( IWS, MINWS );
               if ( LWORK < MINWS ) {

                  // Not enough workspace to use optimal NB: Reduce NB and
                  // determine the minimum value of NB.

                  NB = LWORK / ( SN+1 );
                  NBMIN = max( 2, ILAENV( INBMIN, 'CGEQRF', ' ', SM, SN, -1, -1 ) );


               }
            }
         }

         // Initialize partial column norms. The first N elements of work
         // store the exact column norms.

         for (J = NFXD + 1; J <= N; J++) { // 20
            RWORK( J ) = SCNRM2( SM, A( NFXD+1, J ), 1 );
            RWORK( N+J ) = RWORK( J );
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

               claqps(M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA, JPVT( J ), TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ), WORK( JB+1 ), N-J+1 );

               J = J + FJB;
               GO TO 30;
            }
         } else {
            J = NFXD + 1;
         }

         // Use unblocked code to factor the last or only block.


         if (J <= MINMN) CALL CLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ), TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ) );

      }

      WORK( 1 ) = CMPLX( LWKOPT );
      return;

      // End of CGEQP3

      }

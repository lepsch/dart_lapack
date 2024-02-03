      SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                INB, INBMIN, IXOVER;
      const              INB = 1, INBMIN = 2, IXOVER = 3 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB, NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DLAQP2, DLAQPS, DORMQR, DSWAP, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DNRM2;
      // EXTERNAL ILAENV, DNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments
*  ====================

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }

      if ( INFO == 0 ) {
         MINMN = MIN( M, N )
         if ( MINMN == 0 ) {
            IWS = 1
            LWKOPT = 1
         } else {
            IWS = 3*N + 1
            NB = ILAENV( INB, 'DGEQRF', ' ', M, N, -1, -1 )
            LWKOPT = 2*N + ( N + 1 )*NB
         }
         WORK( 1 ) = LWKOPT

         if ( ( LWORK.LT.IWS ) && .NOT.LQUERY ) {
            INFO = -8
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGEQP3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Move initial columns up front.

      NFXD = 1
      for (J = 1; J <= N; J++) { // 10
         if ( JPVT( J ) != 0 ) {
            if ( J != NFXD ) {
               dswap(M, A( 1, J ), 1, A( 1, NFXD ), 1 );
               JPVT( J ) = JPVT( NFXD )
               JPVT( NFXD ) = J
            } else {
               JPVT( J ) = J
            }
            NFXD = NFXD + 1
         } else {
            JPVT( J ) = J
         }
      } // 10
      NFXD = NFXD - 1

      // Factorize fixed columns
*  =======================

      // Compute the QR factorization of fixed columns and update
      // remaining columns.

      if ( NFXD.GT.0 ) {
         NA = MIN( M, NFXD )
*CC      CALL DGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         dgeqrf(M, NA, A, LDA, TAU, WORK, LWORK, INFO );
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
         if ( NA.LT.N ) {
*CC         CALL DORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
*CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            dormqr('Left', 'Transpose', M, N-NA, NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, LWORK, INFO );
            IWS = MAX( IWS, INT( WORK( 1 ) ) )
         }
      }

      // Factorize free columns
*  ======================

      if ( NFXD.LT.MINMN ) {

         SM = M - NFXD
         SN = N - NFXD
         SMINMN = MINMN - NFXD

         // Determine the block size.

         NB = ILAENV( INB, 'DGEQRF', ' ', SM, SN, -1, -1 )
         NBMIN = 2
         NX = 0

         if ( ( NB.GT.1 ) && ( NB.LT.SMINMN ) ) {

            // Determine when to cross over from blocked to unblocked code.

            NX = MAX( 0, ILAENV( IXOVER, 'DGEQRF', ' ', SM, SN, -1, -1 ) )


            if ( NX.LT.SMINMN ) {

               // Determine if workspace is large enough for blocked code.

               MINWS = 2*SN + ( SN+1 )*NB
               IWS = MAX( IWS, MINWS )
               if ( LWORK.LT.MINWS ) {

                  // Not enough workspace to use optimal NB: Reduce NB and
                  // determine the minimum value of NB.

                  NB = ( LWORK-2*SN ) / ( SN+1 )
                  NBMIN = MAX( 2, ILAENV( INBMIN, 'DGEQRF', ' ', SM, SN, -1, -1 ) )


               }
            }
         }

         // Initialize partial column norms. The first N elements of work
         // store the exact column norms.

         for (J = NFXD + 1; J <= N; J++) { // 20
            WORK( J ) = DNRM2( SM, A( NFXD+1, J ), 1 )
            WORK( N+J ) = WORK( J )
         } // 20

         if ( ( NB.GE.NBMIN ) && ( NB.LT.SMINMN ) && ( NX.LT.SMINMN ) ) {

            // Use blocked code initially.

            J = NFXD + 1

            // Compute factorization: while loop.


            TOPBMN = MINMN - NX
            } // 30
            if ( J.LE.TOPBMN ) {
               JB = MIN( NB, TOPBMN-J+1 )

               // Factorize JB columns among columns J:N.

               dlaqps(M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA, JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ), WORK( 2*N+1 ), WORK( 2*N+JB+1 ), N-J+1 );

               J = J + FJB
               GO TO 30
            }
         } else {
            J = NFXD + 1
         }

         // Use unblocked code to factor the last or only block.


         if (J.LE.MINMN) CALL DLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ), WORK( 2*N+1 ) );

      }

      WORK( 1 ) = IWS
      RETURN

      // End of DGEQP3

      }

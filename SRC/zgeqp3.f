      SUBROUTINE ZGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
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
      // EXTERNAL XERBLA, ZGEQRF, ZLAQP2, ZLAQPS, ZSWAP, ZUNMQR
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DZNRM2;
      // EXTERNAL ILAENV, DZNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments
*  ====================

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }

      if ( INFO.EQ.0 ) {
         MINMN = MIN( M, N )
         if ( MINMN.EQ.0 ) {
            IWS = 1
            LWKOPT = 1
         } else {
            IWS = N + 1
            NB = ILAENV( INB, 'ZGEQRF', ' ', M, N, -1, -1 )
            LWKOPT = ( N + 1 )*NB
         }
         WORK( 1 ) = DCMPLX( LWKOPT )

         if ( ( LWORK.LT.IWS ) .AND. .NOT.LQUERY ) {
            INFO = -8
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZGEQP3', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Move initial columns up front.

      NFXD = 1
      for (J = 1; J <= N; J++) { // 10
         if ( JPVT( J ).NE.0 ) {
            if ( J.NE.NFXD ) {
               zswap(M, A( 1, J ), 1, A( 1, NFXD ), 1 );
               JPVT( J ) = JPVT( NFXD )
               JPVT( NFXD ) = J
            } else {
               JPVT( J ) = J
            }
            NFXD = NFXD + 1
         } else {
            JPVT( J ) = J
         }
   10 CONTINUE
      NFXD = NFXD - 1

      // Factorize fixed columns
*  =======================

      // Compute the QR factorization of fixed columns and update
      // remaining columns.

      if ( NFXD.GT.0 ) {
         NA = MIN( M, NFXD )
*CC      CALL ZGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         zgeqrf(M, NA, A, LDA, TAU, WORK, LWORK, INFO );
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
         if ( NA.LT.N ) {
*CC         CALL ZUNM2R( 'Left', 'Conjugate Transpose', M, N-NA,
*CC  $                   NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK,
*CC  $                   INFO )
            zunmqr('Left', 'Conjugate Transpose', M, N-NA, NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, LWORK, INFO );
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

         NB = ILAENV( INB, 'ZGEQRF', ' ', SM, SN, -1, -1 )
         NBMIN = 2
         NX = 0

         if ( ( NB.GT.1 ) .AND. ( NB.LT.SMINMN ) ) {

            // Determine when to cross over from blocked to unblocked code.

            NX = MAX( 0, ILAENV( IXOVER, 'ZGEQRF', ' ', SM, SN, -1, -1 ) )


            if ( NX.LT.SMINMN ) {

               // Determine if workspace is large enough for blocked code.

               MINWS = ( SN+1 )*NB
               IWS = MAX( IWS, MINWS )
               if ( LWORK.LT.MINWS ) {

                  // Not enough workspace to use optimal NB: Reduce NB and
                  // determine the minimum value of NB.

                  NB = LWORK / ( SN+1 )
                  NBMIN = MAX( 2, ILAENV( INBMIN, 'ZGEQRF', ' ', SM, SN, -1, -1 ) )


               }
            }
         }

         // Initialize partial column norms. The first N elements of work
         // store the exact column norms.

         DO 20 J = NFXD + 1, N
            RWORK( J ) = DZNRM2( SM, A( NFXD+1, J ), 1 )
            RWORK( N+J ) = RWORK( J )
   20    CONTINUE

         if ( ( NB.GE.NBMIN ) .AND. ( NB.LT.SMINMN ) .AND. ( NX.LT.SMINMN ) ) {

            // Use blocked code initially.

            J = NFXD + 1

            // Compute factorization: while loop.


            TOPBMN = MINMN - NX
   30       CONTINUE
            if ( J.LE.TOPBMN ) {
               JB = MIN( NB, TOPBMN-J+1 )

               // Factorize JB columns among columns J:N.

               zlaqps(M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA, JPVT( J ), TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ), WORK( JB+1 ), N-J+1 );

               J = J + FJB
               GO TO 30
            }
         } else {
            J = NFXD + 1
         }

         // Use unblocked code to factor the last or only block.


         IF( J.LE.MINMN ) CALL ZLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ), TAU( J ), RWORK( J ), RWORK( N+J ), WORK( 1 ) )

      }

      WORK( 1 ) = DCMPLX( LWKOPT )
      RETURN

      // End of ZGEQP3

      }

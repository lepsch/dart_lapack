      SUBROUTINE CUNMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NBMAX, LDT, TSIZE;
      const              NBMAX = 64, LDT = NBMAX+1, TSIZE = LDT*NBMAX ;
      // ..
      // .. Local Scalars ..
      bool               LEFT, LQUERY, NOTRAN;
      String             TRANST;
      int                I, I1, I2, I3, IB, IINFO, IWT, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARFB, CLARFT, CUNMR2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK == -1 )

      // NQ is the order of Q and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M
         NW = MAX( 1, N )
      } else {
         NQ = N
         NW = MAX( 1, M )
      }
      if ( .NOT.LEFT && .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 .OR. K.GT.NQ ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, K ) ) {
         INFO = -7
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -10
      } else if ( LWORK.LT.NW && .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO == 0 ) {

         // Compute the workspace requirements

         if ( M == 0 .OR. N == 0 ) {
            LWKOPT = 1
         } else {
            NB = MIN( NBMAX, ILAENV( 1, 'CUNMRQ', SIDE // TRANS, M, N, K, -1 ) )
            LWKOPT = NW*NB + TSIZE
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO != 0 ) {
         xerbla('CUNMRQ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 .OR. N == 0 ) {
         RETURN
      }

      NBMIN = 2
      LDWORK = NW
      if ( NB.GT.1 && NB.LT.K ) {
         if ( LWORK.LT.LWKOPT ) {
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'CUNMRQ', SIDE // TRANS, M, N, K, -1 ) )
         }
      }

      if ( NB.LT.NBMIN .OR. NB.GE.K ) {

         // Use unblocked code

         cunmr2(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO );
      } else {

         // Use blocked code

         IWT = 1 + NW*NB
         if ( ( LEFT && .NOT.NOTRAN ) .OR. ( .NOT.LEFT && NOTRAN ) ) {
            I1 = 1
            I2 = K
            I3 = NB
         } else {
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         }

         if ( LEFT ) {
            NI = N
         } else {
            MI = M
         }

         if ( NOTRAN ) {
            TRANST = 'C'
         } else {
            TRANST = 'N'
         }

         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )

            // Form the triangular factor of the block reflector
            // H = H(i+ib-1) . . . H(i+1) H(i)

            clarft('Backward', 'Rowwise', NQ-K+I+IB-1, IB, A( I, 1 ), LDA, TAU( I ), WORK( IWT ), LDT );
            if ( LEFT ) {

               // H or H**H is applied to C(1:m-k+i+ib-1,1:n)

               MI = M - K + I + IB - 1
            } else {

               // H or H**H is applied to C(1:m,1:n-k+i+ib-1)

               NI = N - K + I + IB - 1
            }

            // Apply H or H**H

            clarfb(SIDE, TRANST, 'Backward', 'Rowwise', MI, NI, IB, A( I, 1 ), LDA, WORK( IWT ), LDT, C, LDC, WORK, LDWORK );
         } // 10
      }
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CUNMRQ

      }

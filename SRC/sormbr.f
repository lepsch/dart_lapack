      SUBROUTINE SORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS, VECT;
      int                INFO, K, LDA, LDC, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               APPLYQ, LEFT, LQUERY, NOTRAN;
      String             TRANST;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORMLQ, SORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )

      // NQ is the order of Q or P and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M
         NW = MAX( 1, N )
      } else {
         NQ = N
         NW = MAX( 1, M )
      }
      if ( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) {
         INFO = -1
      } else if ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -2
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( K.LT.0 ) {
         INFO = -6
      } else if ( ( APPLYQ .AND. LDA.LT.MAX( 1, NQ ) ) .OR. ( .NOT.APPLYQ .AND. LDA.LT.MAX( 1, MIN( NQ, K ) ) ) ) {
         INFO = -8
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -11
      } else if ( LWORK.LT.NW .AND. .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO.EQ.0 ) {
         if ( APPLYQ ) {
            if ( LEFT ) {
               NB = ILAENV( 1, 'SORMQR', SIDE // TRANS, M-1, N, M-1, -1 )
            } else {
               NB = ILAENV( 1, 'SORMQR', SIDE // TRANS, M, N-1, N-1, -1 )
            }
         } else {
            if ( LEFT ) {
               NB = ILAENV( 1, 'SORMLQ', SIDE // TRANS, M-1, N, M-1, -1 )
            } else {
               NB = ILAENV( 1, 'SORMLQ', SIDE // TRANS, M, N-1, N-1, -1 )
            }
         }
         LWKOPT = NW*NB
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         xerbla('SORMBR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      WORK( 1 ) = 1
      if (M.EQ.0 .OR. N.EQ.0) RETURN;

      if ( APPLYQ ) {

         // Apply Q

         if ( NQ.GE.K ) {

            // Q was determined by a call to SGEBRD with nq >= k

            sormqr(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO );
         } else if ( NQ.GT.1 ) {

            // Q was determined by a call to SGEBRD with nq < k

            if ( LEFT ) {
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            } else {
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            }
            sormqr(SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO );
         }
      } else {

         // Apply P

         if ( NOTRAN ) {
            TRANST = 'T'
         } else {
            TRANST = 'N'
         }
         if ( NQ.GT.K ) {

            // P was determined by a call to SGEBRD with nq > k

            sormlq(SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO );
         } else if ( NQ.GT.1 ) {

            // P was determined by a call to SGEBRD with nq <= k

            if ( LEFT ) {
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            } else {
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            }
            sormlq(SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO );
         }
      }
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of SORMBR

      }

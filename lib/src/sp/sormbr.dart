      void sormbr(VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS, VECT;
      int                INFO, K, LDA, LDC, LWORK, M, N;
      double               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               APPLYQ, LEFT, LQUERY, NOTRAN;
      String             TRANST;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORMLQ, SORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      APPLYQ = lsame( VECT, 'Q' );
      LEFT = lsame( SIDE, 'L' );
      NOTRAN = lsame( TRANS, 'N' );
      LQUERY = ( LWORK == -1 );

      // NQ is the order of Q or P and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M;
         NW = max( 1, N );
      } else {
         NQ = N;
         NW = max( 1, M );
      }
      if ( !APPLYQ && !lsame( VECT, 'P' ) ) {
         INFO = -1;
      } else if ( !LEFT && !lsame( SIDE, 'R' ) ) {
         INFO = -2;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( K < 0 ) {
         INFO = -6;
      } else if ( ( APPLYQ && LDA < max( 1, NQ ) ) || ( !APPLYQ && LDA < max( 1, min( NQ, K ) ) ) ) {
         INFO = -8;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -11;
      } else if ( LWORK < NW && !LQUERY ) {
         INFO = -13;
      }

      if ( INFO == 0 ) {
         if ( APPLYQ ) {
            if ( LEFT ) {
               NB = ilaenv( 1, 'SORMQR', SIDE + TRANS, M-1, N, M-1, -1 );
            } else {
               NB = ilaenv( 1, 'SORMQR', SIDE + TRANS, M, N-1, N-1, -1 );
            }
         } else {
            if ( LEFT ) {
               NB = ilaenv( 1, 'SORMLQ', SIDE + TRANS, M-1, N, M-1, -1 );
            } else {
               NB = ilaenv( 1, 'SORMLQ', SIDE + TRANS, M, N-1, N-1, -1 );
            }
         }
         LWKOPT = NW*NB;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('SORMBR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      WORK[1] = 1;
      if (M == 0 || N == 0) return;

      if ( APPLYQ ) {

         // Apply Q

         if ( NQ >= K ) {

            // Q was determined by a call to SGEBRD with nq >= k

            sormqr(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO );
         } else if ( NQ > 1 ) {

            // Q was determined by a call to SGEBRD with nq < k

            if ( LEFT ) {
               MI = M - 1;
               NI = N;
               I1 = 2;
               I2 = 1;
            } else {
               MI = M;
               NI = N - 1;
               I1 = 1;
               I2 = 2;
            }
            sormqr(SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO );
         }
      } else {

         // Apply P

         if ( NOTRAN ) {
            TRANST = 'T';
         } else {
            TRANST = 'N';
         }
         if ( NQ > K ) {

            // P was determined by a call to SGEBRD with nq > k

            sormlq(SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO );
         } else if ( NQ > 1 ) {

            // P was determined by a call to SGEBRD with nq <= k

            if ( LEFT ) {
               MI = M - 1;
               NI = N;
               I1 = 2;
               I2 = 1;
            } else {
               MI = M;
               NI = N - 1;
               I1 = 1;
               I2 = 2;
            }
            sormlq(SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO );
         }
      }
      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

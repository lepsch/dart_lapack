      void cunmqr(SIDE, TRANS, M, N, K, final Matrix<double> A, final int LDA, TAU, final Matrix<double> C, final int LDC, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, LWORK, M, N;
      Complex            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

      int                NBMAX, LDT, TSIZE;
      const              NBMAX = 64, LDT = NBMAX+1, TSIZE = LDT*NBMAX ;
      bool               LEFT, LQUERY, NOTRAN;
      int                I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARFB, CLARFT, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      LEFT = lsame( SIDE, 'L' );
      NOTRAN = lsame( TRANS, 'N' );
      LQUERY = ( LWORK == -1 );

      // NQ is the order of Q and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M;
         NW = max( 1, N );
      } else {
         NQ = N;
         NW = max( 1, M );
      }
      if ( !LEFT && !lsame( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 || K > NQ ) {
         INFO = -5;
      } else if ( LDA < max( 1, NQ ) ) {
         INFO = -7;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -10;
      } else if ( LWORK < NW && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO == 0 ) {

         // Compute the workspace requirements

         NB = min( NBMAX, ilaenv( 1, 'CUNMQR', SIDE + TRANS, M, N, K, -1 ) );
         LWKOPT = NW*NB + TSIZE;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('CUNMQR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 || K == 0 ) {
         WORK[1] = 1;
         return;
      }

      NBMIN = 2;
      LDWORK = NW;
      if ( NB > 1 && NB < K ) {
         if ( LWORK < LWKOPT ) {
            NB = (LWORK-TSIZE) / LDWORK;
            NBMIN = max( 2, ilaenv( 2, 'CUNMQR', SIDE + TRANS, M, N, K, -1 ) );
         }
      }

      if ( NB < NBMIN || NB >= K ) {

         // Use unblocked code

         cunm2r(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO );
      } else {

         // Use blocked code

         IWT = 1 + NW*NB;
         if ( ( LEFT && !NOTRAN ) || ( !LEFT && NOTRAN ) ) {
            I1 = 1;
            I2 = K;
            I3 = NB;
         } else {
            I1 = ( ( K-1 ) / NB )*NB + 1;
            I2 = 1;
            I3 = -NB;
         }

         if ( LEFT ) {
            NI = N;
            JC = 1;
         } else {
            MI = M;
            IC = 1;
         }

         for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) { // 10
            IB = min( NB, K-I+1 );

            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1)

            clarft('Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), LDA, TAU( I ), WORK( IWT ), LDT );
            if ( LEFT ) {

               // H or H**H is applied to C(i:m,1:n)

               MI = M - I + 1;
               IC = I;
            } else {

               // H or H**H is applied to C(1:m,i:n)

               NI = N - I + 1;
               JC = I;
            }

            // Apply H or H**H

            clarfb(SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, IB, A( I, I ), LDA, WORK( IWT ), LDT, C( IC, JC ), LDC, WORK, LDWORK );
         } // 10
      }
      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

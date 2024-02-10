      void sormtr(final int SIDE, final int UPLO, final int TRANS, final int M, final int N, final Matrix<double> A, final int LDA, final int TAU, final Matrix<double> C, final int LDC, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS, UPLO;
      int                INFO, LDA, LDC, LWORK, M, N;
      double               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LEFT, LQUERY, UPPER;
      int                I1, I2, IINFO, LWKOPT, MI, NI, NB, NQ, NW;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORMQL, SORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input arguments

      INFO = 0;
      LEFT = lsame( SIDE, 'L' );
      UPPER = lsame( UPLO, 'U' );
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
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( !lsame( TRANS, 'N' ) && !lsame( TRANS, 'T' ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, NQ ) ) {
         INFO = -7;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -10;
      } else if ( LWORK < NW && !LQUERY ) {
         INFO = -12;
      }

      if ( INFO == 0 ) {
         if ( UPPER ) {
            if ( LEFT ) {
               NB = ilaenv( 1, 'SORMQL', SIDE + TRANS, M-1, N, M-1, -1 );
            } else {
               NB = ilaenv( 1, 'SORMQL', SIDE + TRANS, M, N-1, N-1, -1 );
            }
         } else {
            if ( LEFT ) {
               NB = ilaenv( 1, 'SORMQR', SIDE + TRANS, M-1, N, M-1, -1 );
            } else {
               NB = ilaenv( 1, 'SORMQR', SIDE + TRANS, M, N-1, N-1, -1 );
            }
         }
         LWKOPT = NW*NB;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('SORMTR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 || NQ == 1 ) {
         WORK[1] = 1;
         return;
      }

      if ( LEFT ) {
         MI = M - 1;
         NI = N;
      } else {
         MI = M;
         NI = N - 1;
      }

      if ( UPPER ) {

         // Q was determined by a call to SSYTRD with UPLO = 'U'

         sormql(SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, LDC, WORK, LWORK, IINFO );
      } else {

         // Q was determined by a call to SSYTRD with UPLO = 'L'

         if ( LEFT ) {
            I1 = 2;
            I2 = 1;
         } else {
            I1 = 1;
            I2 = 2;
         }
         sormqr(SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO );
      }
      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

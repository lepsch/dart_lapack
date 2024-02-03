      SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NBMAX, LDT, TSIZE;
      const              NBMAX = 64, LDT = NBMAX+1, TSIZE = LDT*NBMAX ;
      // ..
      // .. Local Scalars ..
      bool               LEFT, LQUERY, NOTRAN;
      int                I, I1, I2, I3, IB, IINFO, IWT, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB, DLARFT, DORM2L, XERBLA
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
      } else if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -2
      } else if ( M < 0 ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( K < 0 || K > NQ ) {
         INFO = -5
      } else if ( LDA < MAX( 1, NQ ) ) {
         INFO = -7
      } else if ( LDC < MAX( 1, M ) ) {
         INFO = -10
      } else if ( LWORK < NW && .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO == 0 ) {

         // Compute the workspace requirements

         if ( M == 0 || N == 0 ) {
            LWKOPT = 1
         } else {
            NB = MIN( NBMAX, ILAENV( 1, 'DORMQL', SIDE // TRANS, M, N, K, -1 ) )
            LWKOPT = NW*NB + TSIZE
         }
         WORK( 1 ) = LWKOPT
      }

      if ( INFO != 0 ) {
         xerbla('DORMQL', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         RETURN
      }

      NBMIN = 2
      LDWORK = NW
      if ( NB > 1 && NB < K ) {
         if ( LWORK < LWKOPT ) {
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQL', SIDE // TRANS, M, N, K, -1 ) )
         }
      }

      if ( NB < NBMIN || NB >= K ) {

         // Use unblocked code

         dorm2l(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO );
      } else {

         // Use blocked code

         IWT = 1 + NW*NB
         if ( ( LEFT && NOTRAN ) || ( .NOT.LEFT && .NOT.NOTRAN ) ) {
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

         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )

            // Form the triangular factor of the block reflector
            // H = H(i+ib-1) . . . H(i+1) H(i)

            dlarft('Backward', 'Columnwise', NQ-K+I+IB-1, IB, A( 1, I ), LDA, TAU( I ), WORK( IWT ), LDT );
            if ( LEFT ) {

               // H or H**T is applied to C(1:m-k+i+ib-1,1:n)

               MI = M - K + I + IB - 1
            } else {

               // H or H**T is applied to C(1:m,1:n-k+i+ib-1)

               NI = N - K + I + IB - 1
            }

            // Apply H or H**T

            dlarfb(SIDE, TRANS, 'Backward', 'Columnwise', MI, NI, IB, A( 1, I ), LDA, WORK( IWT ), LDT, C, LDC, WORK, LDWORK );
         } // 10
      }
      WORK( 1 ) = LWKOPT
      RETURN

      // End of DORMQL

      }

      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

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
      int                I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFB, DLARFT, DORM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )

      // NQ is the order of Q and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M
         NW = MAX( 1, N )
      } else {
         NQ = N
         NW = MAX( 1, M )
      }
      if ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 .OR. K.GT.NQ ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, NQ ) ) {
         INFO = -7
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -10
      } else if ( LWORK.LT.NW .AND. .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO.EQ.0 ) {

         // Compute the workspace requirements

         NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K, -1 ) )
         LWKOPT = NW*NB + TSIZE
         WORK( 1 ) = LWKOPT
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      NBMIN = 2
      LDWORK = NW
      if ( NB.GT.1 .AND. NB.LT.K ) {
         if ( LWORK.LT.LWKOPT ) {
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K, -1 ) )
         }
      }

      if ( NB.LT.NBMIN .OR. NB.GE.K ) {

         // Use unblocked code

         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )
      } else {

         // Use blocked code

         IWT = 1 + NW*NB
         if ( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) {
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
            JC = 1
         } else {
            MI = M
            IC = 1
         }

         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )

            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1)

            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), LDA, TAU( I ), WORK( IWT ), LDT )
            if ( LEFT ) {

               // H or H**T is applied to C(i:m,1:n)

               MI = M - I + 1
               IC = I
            } else {

               // H or H**T is applied to C(1:m,i:n)

               NI = N - I + 1
               JC = I
            }

            // Apply H or H**T

            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, IB, A( I, I ), LDA, WORK( IWT ), LDT, C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      }
      WORK( 1 ) = LWKOPT
      RETURN

      // End of DORMQR

      }

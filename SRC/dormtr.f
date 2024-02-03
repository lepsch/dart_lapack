      SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS, UPLO;
      int                INFO, LDA, LDC, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LEFT, LQUERY, UPPER;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DORMQL, DORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
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
      } else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, NQ ) ) {
         INFO = -7
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -10
      } else if ( LWORK.LT.NW .AND. .NOT.LQUERY ) {
         INFO = -12
      }

      if ( INFO.EQ.0 ) {
         if ( UPPER ) {
            if ( LEFT ) {
               NB = ILAENV( 1, 'DORMQL', SIDE // TRANS, M-1, N, M-1, -1 )
            } else {
               NB = ILAENV( 1, 'DORMQL', SIDE // TRANS, M, N-1, N-1, -1 )
            }
         } else {
            if ( LEFT ) {
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1, -1 )
            } else {
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1, -1 )
            }
         }
         LWKOPT = NW*NB
         WORK( 1 ) = LWKOPT
      }

      if ( INFO.NE.0 ) {
         xerbla('DORMTR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( LEFT ) {
         MI = M - 1
         NI = N
      } else {
         MI = M
         NI = N - 1
      }

      if ( UPPER ) {

         // Q was determined by a call to DSYTRD with UPLO = 'U'

         dormql(SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, LDC, WORK, LWORK, IINFO );
      } else {

         // Q was determined by a call to DSYTRD with UPLO = 'L'

         if ( LEFT ) {
            I1 = 2
            I2 = 1
         } else {
            I1 = 1
            I2 = 2
         }
         dormqr(SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO );
      }
      WORK( 1 ) = LWKOPT
      RETURN

      // End of DORMTR

      }

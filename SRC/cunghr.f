      SUBROUTINE CUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LWKOPT, NB, NH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CUNGQR, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      NH = IHI - ILO
      LQUERY = ( LWORK.EQ.-1 )
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) {
         INFO = -2
      } else if ( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LWORK.LT.MAX( 1, NH ) .AND. .NOT.LQUERY ) {
         INFO = -8
      }

      if ( INFO.EQ.0 ) {
         NB = ILAENV( 1, 'CUNGQR', ' ', NH, NH, NH, -1 )
         LWKOPT = MAX( 1, NH )*NB
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         xerbla('CUNGHR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      // Shift the vectors which define the elementary reflectors one
      // column to the right, and set the first ilo and the last n-ihi
      // rows and columns to those of the unit matrix

      DO 40 J = IHI, ILO + 1, -1
         for (I = 1; I <= J - 1; I++) { // 10
            A( I, J ) = ZERO
         } // 10
         for (I = J + 1; I <= IHI; I++) { // 20
            A( I, J ) = A( I, J-1 )
         } // 20
         for (I = IHI + 1; I <= N; I++) { // 30
            A( I, J ) = ZERO
         } // 30
      } // 40
      for (J = 1; J <= ILO; J++) { // 60
         for (I = 1; I <= N; I++) { // 50
            A( I, J ) = ZERO
         } // 50
         A( J, J ) = ONE
      } // 60
      for (J = IHI + 1; J <= N; J++) { // 80
         for (I = 1; I <= N; I++) { // 70
            A( I, J ) = ZERO
         } // 70
         A( J, J ) = ONE
      } // 80

      if ( NH.GT.0 ) {

         // Generate Q(ilo+1:ihi,ilo+1:ihi)

         cungqr(NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), WORK, LWORK, IINFO );
      }
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CUNGHR

      }

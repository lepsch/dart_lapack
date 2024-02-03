      SUBROUTINE CUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
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
      bool               LQUERY, UPPER;
      int                I, IINFO, J, LWKOPT, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CUNGQL, CUNGQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) {
         INFO = -7
      }

      if ( INFO.EQ.0 ) {
         if ( UPPER ) {
           NB = ILAENV( 1, 'CUNGQL', ' ', N-1, N-1, N-1, -1 )
         } else {
           NB = ILAENV( 1, 'CUNGQR', ' ', N-1, N-1, N-1, -1 )
         }
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         xerbla('CUNGTR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( UPPER ) {

         // Q was determined by a call to CHETRD with UPLO = 'U'

         // Shift the vectors which define the elementary reflectors one
         // column to the left, and set the last row and column of Q to
         // those of the unit matrix

         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
            } // 10
            A( N, J ) = ZERO
         } // 20
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
         } // 30
         A( N, N ) = ONE

         // Generate Q(1:n-1,1:n-1)

         cungql(N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO );

      } else {

         // Q was determined by a call to CHETRD with UPLO = 'L'.

         // Shift the vectors which define the elementary reflectors one
         // column to the right, and set the first row and column of Q to
         // those of the unit matrix

         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
            } // 40
         } // 50
         A( 1, 1 ) = ONE
         for (I = 2; I <= N; I++) { // 60
            A( I, 1 ) = ZERO
         } // 60
         if ( N.GT.1 ) {

            // Generate Q(2:n,2:n)

            cungqr(N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO );
         }
      }
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CUNGTR

      }

      SUBROUTINE CUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             VECT;
      int                INFO, K, LDA, LWORK, M, N;
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
      bool               LQUERY, WANTQ;
      int                I, IINFO, J, LWKOPT, MN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CUNGLQ, CUNGQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M, K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT. MIN( N, K ) ) ) ) {
         INFO = -3
      } else if ( K.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -6
      } else if ( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) {
         INFO = -9
      }

      if ( INFO.EQ.0 ) {
         WORK( 1 ) = 1
         if ( WANTQ ) {
            if ( M.GE.K ) {
               CALL CUNGQR( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
            } else {
               if ( M.GT.1 ) {
                  CALL CUNGQR( M-1, M-1, M-1, A, LDA, TAU, WORK, -1, IINFO )
               }
            }
         } else {
            if ( K.LT.N ) {
               CALL CUNGLQ( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
            } else {
               if ( N.GT.1 ) {
                  CALL CUNGLQ( N-1, N-1, N-1, A, LDA, TAU, WORK, -1, IINFO )
               }
            }
         }
         LWKOPT = INT( WORK( 1 ) )
         LWKOPT = MAX (LWKOPT, MN)
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CUNGBR', -INFO )
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( WANTQ ) {

         // Form Q, determined by a call to CGEBRD to reduce an m-by-k
         // matrix

         if ( M.GE.K ) {

            // If m >= k, assume m >= n >= k

            CALL CUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )

         } else {

            // If m < k, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // column to the right, and set the first row and column of Q
           t // o those of the unit matrix

            DO 20 J = M, 2, -1
               A( 1, J ) = ZERO
               DO 10 I = J + 1, M
                  A( I, J ) = A( I, J-1 )
   10          CONTINUE
   20       CONTINUE
            A( 1, 1 ) = ONE
            DO 30 I = 2, M
               A( I, 1 ) = ZERO
   30       CONTINUE
            if ( M.GT.1 ) {

               // Form Q(2:m,2:m)

               CALL CUNGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
            }
         }
      } else {

         // Form P**H, determined by a call to CGEBRD to reduce a k-by-n
         // matrix

         if ( K.LT.N ) {

            // If k < n, assume k <= m <= n

            CALL CUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )

         } else {

            // If k >= n, assume m = n

            // Shift the vectors which define the elementary reflectors one
            // row downward, and set the first row and column of P**H to
           t // hose of the unit matrix

            A( 1, 1 ) = ONE
            DO 40 I = 2, N
               A( I, 1 ) = ZERO
   40       CONTINUE
            DO 60 J = 2, N
               DO 50 I = J - 1, 2, -1
                  A( I, J ) = A( I-1, J )
   50          CONTINUE
               A( 1, J ) = ZERO
   60       CONTINUE
            if ( N.GT.1 ) {

               // Form P**H(2:n,2:n)

               CALL CUNGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
            }
         }
      }
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CUNGBR

      }

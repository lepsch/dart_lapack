      SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          VECT
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WANTQ
      INTEGER            I, IINFO, J, LWKOPT, MN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORGLQ, DORGQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M, K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT. MIN( N, K ) ) ) ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = 1
         IF( WANTQ ) THEN
            IF( M.GE.K ) THEN
               CALL DORGQR( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
            ELSE
               IF( M.GT.1 ) THEN
                  CALL DORGQR( M-1, M-1, M-1, A, LDA, TAU, WORK, -1, IINFO )
               END IF
            END IF
         ELSE
            IF( K.LT.N ) THEN
               CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
            ELSE
               IF( N.GT.1 ) THEN
                  CALL DORGLQ( N-1, N-1, N-1, A, LDA, TAU, WORK, -1, IINFO )
               END IF
            END IF
         END IF
         LWKOPT = INT( WORK( 1 ) )
         LWKOPT = MAX (LWKOPT, MN)
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK( 1 ) = LWKOPT
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( WANTQ ) THEN
*
*        Form Q, determined by a call to DGEBRD to reduce an m-by-k
*        matrix
*
         IF( M.GE.K ) THEN
*
*           If m >= k, assume m >= n >= k
*
            CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*           If m < k, assume m = n
*
*           Shift the vectors which define the elementary reflectors one
*           column to the right, and set the first row and column of Q
*           to those of the unit matrix
*
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
            IF( M.GT.1 ) THEN
*
*              Form Q(2:m,2:m)
*
               CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
            END IF
         END IF
      ELSE
*
*        Form P**T, determined by a call to DGEBRD to reduce a k-by-n
*        matrix
*
         IF( K.LT.N ) THEN
*
*           If k < n, assume k <= m <= n
*
            CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*           If k >= n, assume m = n
*
*           Shift the vectors which define the elementary reflectors one
*           row downward, and set the first row and column of P**T to
*           those of the unit matrix
*
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
            IF( N.GT.1 ) THEN
*
*              Form P**T(2:n,2:n)
*
               CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
            END IF
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORGBR
*
      END

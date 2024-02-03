      SUBROUTINE SORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LWKOPT, NB, NH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SORGQR, XERBLA
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
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, NH ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF

      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'SORGQR', ' ', NH, NH, NH, -1 )
         LWKOPT = MAX( 1, NH )*NB
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF

      // Shift the vectors which define the elementary reflectors one
      // column to the right, and set the first ilo and the last n-ihi
      // rows and columns to those of the unit matrix

      DO 40 J = IHI, ILO + 1, -1
         DO 10 I = 1, J - 1
            A( I, J ) = ZERO
   10    CONTINUE
         DO 20 I = J + 1, IHI
            A( I, J ) = A( I, J-1 )
   20    CONTINUE
         DO 30 I = IHI + 1, N
            A( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1, ILO
         DO 50 I = 1, N
            A( I, J ) = ZERO
   50    CONTINUE
         A( J, J ) = ONE
   60 CONTINUE
      DO 80 J = IHI + 1, N
         DO 70 I = 1, N
            A( I, J ) = ZERO
   70    CONTINUE
         A( J, J ) = ONE
   80 CONTINUE

      IF( NH.GT.0 ) THEN

         // Generate Q(ilo+1:ihi,ilo+1:ihi)

         CALL SORGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of SORGHR

      }
